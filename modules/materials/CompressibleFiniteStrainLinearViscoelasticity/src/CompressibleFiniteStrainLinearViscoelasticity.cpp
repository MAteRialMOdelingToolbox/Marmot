#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  double CompressibleFiniteStrainLinearViscoelasticity::getDensity( const double* stateVars ) const
  {

    if ( nMaterialProperties <
         ( 2 + nElasticPropertiesMap.at( hyperelasticBase ) + 1 + maxwellProperties.nMaxwell * 2 ) ) {
      throw std::runtime_error(
        "CompressibleFiniteStrainLinearViscoelasticity::getDensity: Not enough material properties provided." );
    }
    else
      return materialProperties[nMaterialProperties - 1];
  }

  std::tuple< double,
              FastorStandardTensors::Tensor33d,
              FastorStandardTensors::Tensor3333d,
              FastorStandardTensors::Tensor333333d >
  CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives(
    const FastorStandardTensors::Tensor33d& C ) const
  {

    // select hyperelastic base
    // compute enregy density and derivatives
    // with autodiff
    // we need a function mapping from tensor33 with dual3rd to dual3rd
    std::function< autodiff::dual3rd( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& ) >
      energyDensityFunction;

    switch ( hyperelasticBase ) {
    case NeoHooke:
      return ContinuumMechanics::EnergyDensityFunctions::Compressible::ThirdOrderDerived::
        standardNeoHooke( C, elasticProperties[0], elasticProperties[1] );
    case PenceGouNeoHooke:
      energyDensityFunction = [this]( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& C_ad ) {
        autodiff::dual3rd res = ContinuumMechanics::EnergyDensityFunctions::Compressible::PenceGouPotentialB<
          autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1] );
        return res;
      };
      break;
    case Yeoh:
      energyDensityFunction = [this]( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& C_ad ) {
        autodiff::dual3rd res = ContinuumMechanics::EnergyDensityFunctions::Compressible::YeohPotential<
          autodiff::dual3rd >( C_ad,
                               elasticProperties[0],
                               elasticProperties[1],
                               elasticProperties[2],
                               elasticProperties[3] );
        return res;
      };
      break;
    case MooneyRivlin:
      energyDensityFunction = [this]( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& C_ad ) {
        autodiff::dual3rd res = ContinuumMechanics::EnergyDensityFunctions::Compressible::MooneyRivlinPotential<
          autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1], elasticProperties[2] );
        return res;
      };
      break;
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives: "
                                "Unknown hyperelastic base." );
    }
    return Marmot::AutomaticDifferentiation::ThirdOrder::d3f_dT3< 3 >( energyDensityFunction, C );
  }

  CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      hyperelasticBase( static_cast< HyperelasticBase >( materialProperties[0] ) ),
      onlyShearCreep( materialProperties[1] ),
      elasticProperties( &materialProperties[2], nElasticPropertiesMap.at( hyperelasticBase ) ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::
          createMaxwellProperties( materialProperties[2 + nElasticPropertiesMap.at( hyperelasticBase )],
                                   &materialProperties[3 + nElasticPropertiesMap.at( hyperelasticBase )] ) )
  {

    if ( onlyShearCreep == 1.0 ) {
      double G = 0;
      switch ( hyperelasticBase ) {
      case NeoHooke:
      case PenceGouNeoHooke: G = elasticProperties[1]; break;
      case Yeoh:
      case MooneyRivlin: G = 2 * elasticProperties[0]; break;
      default:
        throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::"
                                  "CompressibleFiniteStrainLinearViscoelasticity: Unknown hyperelastic base." );
      }

      initialCompliance = .5 / G * Spatial3D::I4;
    }

    else {

      // start from initial stiffness tensor dPK2 / dE,
      Tensor3333d initialStiffness = 4. *
                                     std::get< 2 >( computeEnergyDensityAndDerivatives( evaluate( Spatial3D::I ) ) );

      initialCompliance = invertMinorSymmetricFourthOrderTensor( initialStiffness );
    }

    initializeStateLayout();
  }

  void CompressibleFiniteStrainLinearViscoelasticity::computeStress( ConstitutiveResponse< 3 >& response,
                                                                     AlgorithmicModuli< 3 >&    tangents,
                                                                     const Deformation< 3 >&    deformation,
                                                                     const TimeIncrement&       timeIncrement ) const
  {
    // map to deformation gradient
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F );

    // compute energy density, first, second and third partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = computeEnergyDensityAndDerivatives( C );

    // compute initial PK2 stress
    const Tensor33d PK2zero = 2. * dPsi_dC;

    // split into volumetric and deviatoric parts if onlyShearCreep is 1
    // otherwise keep PK2dev representing the full PK2 stress
    const Tensor33d PK2vol = ( onlyShearCreep / 3.0 ) * trace( PK2zero ) * Spatial3D::I;
    Tensor33d       PK2dev = PK2zero - PK2vol;

    // conpute increment in initial PK2dev
    TensorMap33d    PK2dev_old = stateLayout.getAs< TensorMap33d >( response.stateVars, "S0_old" );
    const Tensor33d dPK2dev    = PK2dev - PK2dev_old;

    // update old PK2dev for next increment
    memcpy( PK2dev_old.data(), PK2dev.data(), 9 * sizeof( double ) );

    // tangents needed for viscoelastic contribution
    const Tensor3333d dPK2_dE    = 4. * d2Psi_dCdC;
    const Tensor3333d dPK2vol_dE = onlyShearCreep / 3.0 *
                                   outer( Spatial3D::I, einsum< ijkl, kl >( dPK2_dE, Spatial3D::I ) );
    Tensor3333d         dPK2dev_dE    = dPK2_dE - dPK2vol_dE;
    const Tensor333333d d2PK2_dEdE    = 8. * d3Psi_dCdCdC;
    const Tensor333333d d2PK2vol_dEdE = onlyShearCreep / 3 *
                                        einsum< ij, mnKL >( Spatial3D::I,
                                                            einsum< ijklmn, ij >( d2PK2_dEdE, Spatial3D::I ) );
    const Tensor333333d d2PK2dev_dEdE = d2PK2_dEdE - d2PK2vol_dEdE;

    // add viscoelastic contribution to (deviatoric if onlyShearCreep) PK2 stress
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel( PK2dev,
                                                                                        dPK2dev_dE,
                                                                                        d2PK2dev_dEdE,
                                                                                        initialCompliance,
                                                                                        dPK2dev,
                                                                                        timeIncrement.dT,
                                                                                        maxwellProperties,
                                                                                        stateLayout
                                                                                          .getPtr( response.stateVars,
                                                                                                   "creepStateVars" ) );
    // compute total PK2 stress
    const Tensor33d PK2 = evaluate( PK2vol + PK2dev );

    // compute Kirchhoff stress from PK2 stress
    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F );
    response.tau                         = tau;
    response.elasticEnergyDensity        = psi_; // TODO: replace by actual elastic energy density

    // compute tangent operator dtau / dF
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, evaluate( dPK2vol_dE + dPK2dev_dE ) ),
                                             0.5 * dC_dF ) +
                       dTau_dF;
  }
} // namespace Marmot::Materials
