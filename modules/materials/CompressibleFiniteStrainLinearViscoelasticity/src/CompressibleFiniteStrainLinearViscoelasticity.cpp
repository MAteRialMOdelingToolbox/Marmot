#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <functional>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace Marmot::Math::Tensors::FastorIndices;
  using namespace Marmot::Math::Tensors::FastorTensors;

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
              Marmot::Math::Tensors::FastorTensors::Tensor33d,
              Marmot::Math::Tensors::FastorTensors::Tensor3333d,
              Marmot::Math::Tensors::FastorTensors::Tensor333333d >
  CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives(
    const Marmot::Math::Tensors::FastorTensors::Tensor33d& C ) const
  {

    // select hyperelastic base
    // compute enregy density and derivatives
    // with autodiff
    // we need a function mapping from tensor33 with dual3rd to dual3rd
    std::function< autodiff::dual3rd( const Marmot::Math::Tensors::FastorTensors::Tensor33t< autodiff::dual3rd >& ) >
      energyDensityFunction;

    // every base except PenceGouNeoHooke composes Psi = Psi_iso(shape params) + Psi_vol(kappa),
    // with kappa the LAST elasticProperties entry; PenceGouNeoHooke keeps its own (K, G) potential
    // -- see the class docs.
    namespace EDF               = ContinuumMechanics::Elasticity::EnergyDensityFunctions;
    const auto&  par            = elasticProperties;
    const double kappa          = par( par.size() - 1 );
    const auto   volumetricPart = [kappa]( const auto& C_ ) { return EDF::Volumetric::PenaltyPotential( C_, kappa ); };

    switch ( hyperelasticBase ) {
    case NeoHooke:
      energyDensityFunction = EDF::
        makeDecoupledPotential( [mu = par( 0 )](
                                  const auto& C_ ) { return EDF::Isochoric::NeoHookePotential( C_, mu ); },
                                volumetricPart );
      break;
    case Yeoh:
      energyDensityFunction = EDF::
        makeDecoupledPotential( [c10 = par( 0 ), c20 = par( 1 ), c30 = par( 2 )](
                                  const auto& C_ ) { return EDF::Isochoric::YeohPotential( C_, c10, c20, c30 ); },
                                volumetricPart );
      break;
    case PenceGouNeoHooke:
      energyDensityFunction = [K = par( 0 ), G = par( 1 )]( const Tensor33t< autodiff::dual3rd >& C_ad ) {
        return EDF::PenceGouPotentialB( C_ad, K, G );
      };
      break;
    case MooneyRivlin:
      energyDensityFunction = EDF::
        makeDecoupledPotential( [c10 = par( 0 ), c01 = par( 1 )](
                                  const auto& C_ ) { return EDF::Isochoric::MooneyRivlinPotential( C_, c10, c01 ); },
                                volumetricPart );
      break;
    case ArrudaBoyce:
      energyDensityFunction = EDF::
        makeDecoupledPotential( [mu = par( 0 ), lambdaL = par( 1 )](
                                  const auto&
                                    C_ ) { return EDF::Isochoric::ArrudaBoyce8ChainPotential( C_, mu, lambdaL ); },
                                volumetricPart );
      break;
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives: "
                                "Unknown hyperelastic base." );
    }
    auto [psi,
          dPsi_dC,
          d2Psi_dCdC,
          d3Psi_dCdCdC] = Marmot::Math::Differentiation::Automatic::ThirdOrder::d3f_dT3< 3 >( energyDensityFunction,
                                                                                              C,
                                                                                              true );

    // The potentials are evaluated on a general (not necessarily symmetric) C, so the
    // autodiff derivatives w.r.t. the independent entries C_ij and C_ji differ (e.g.
    // tr(C C) gives d2/dC01dC01 = 0 but d2/dC01dC10 = 2). Only their symmetric part acts on
    // symmetric increments of C, so restore the minor symmetries here; without this, the
    // stiffness inverted for onlyShearCreep = 0 is singular in the shear block.
    Tensor3333d   d2Sym;
    Tensor333333d d3Sym;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ ) {
            d2Sym( i, j, k, l ) = 0.25 * ( d2Psi_dCdC( i, j, k, l ) + d2Psi_dCdC( j, i, k, l ) +
                                           d2Psi_dCdC( i, j, l, k ) + d2Psi_dCdC( j, i, l, k ) );
            for ( int m = 0; m < 3; m++ )
              for ( int n = 0; n < 3; n++ ) {
                double sum = 0.0;
                for ( int a = 0; a < 2; a++ )
                  for ( int b = 0; b < 2; b++ )
                    for ( int c = 0; c < 2; c++ )
                      sum += d3Psi_dCdCdC( a ? j : i, a ? i : j, b ? l : k, b ? k : l, c ? n : m, c ? m : n );
                d3Sym( i, j, k, l, m, n ) = 0.125 * sum;
              }
          }

    return { psi, dPsi_dC, d2Sym, d3Sym };
  }

  CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : Marmot::Materials::MaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      hyperelasticBase( static_cast< HyperelasticBase >( materialProperties[0] ) ),
      onlyShearCreep( materialProperties[1] ),
      elasticProperties( &materialProperties[2], nElasticPropertiesMap.at( hyperelasticBase ) ),
      maxwellProperties(
        ContinuumMechanics::Viscoelasticity::
          createMaxwellProperties( static_cast< int >(
                                     materialProperties[2 + nElasticPropertiesMap.at( hyperelasticBase )] ),
                                   &materialProperties[3 + nElasticPropertiesMap.at( hyperelasticBase )] ) )
  {

    if ( hyperelasticBase == ArrudaBoyce && !( elasticProperties[1] > 1.0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": ArrudaBoyce locking stretch lambdaL must be > 1, got "
                                   << elasticProperties[1] );

    if ( onlyShearCreep == 1.0 ) {
      double G = 0;
      switch ( hyperelasticBase ) {
      case NeoHooke: G = elasticProperties[0]; break;
      case Yeoh: G = 2 * elasticProperties[0]; break;
      // reference shear modulus of the isochoric Mooney-Rivlin potential is 2*(C10+C01), not
      // 2*C10 alone -- C01 also contributes to dPsi/dIbar1 at the reference configuration via
      // the Ibar2 term (dIbar2/dIbar1|C=I effectively contributes an equal weight to C01 there).
      case MooneyRivlin: G = 2 * ( elasticProperties[0] + elasticProperties[1] ); break;
      case PenceGouNeoHooke: G = elasticProperties[1]; break;
      // mu is NOT the infinitesimal shear modulus at finite lambdaL: G0 = 2*dPsi_iso/dIbar1
      // evaluated at Ibar1=3 (the reference configuration), which for the closed-form
      // Arruda-Boyce potential (see Detail::arrudaBoyce8ChainEnergyAndDerivative) works out to
      // mu/3 + 2*mu/(3*(1-1/lambdaL^2)); only in the lambdaL->infinity limit does this reduce to mu.
      case ArrudaBoyce: {
        const double lambdaL = elasticProperties[1];
        const double w0      = 1. - 1. / ( lambdaL * lambdaL );
        G                    = elasticProperties[0] / 3. + 2. * elasticProperties[0] / ( 3. * w0 );
        break;
      }
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

      initialCompliance = Marmot::Math::Tensors::invertMinorSymmetricFourthOrderTensor( initialStiffness );
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
    const auto [C,
                dC_dF] = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::FirstOrder::rightCauchyGreen( F );

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
    ContinuumMechanics::Viscoelasticity::evaluateGeneralizedMaxwellModel( PK2dev,
                                                                          dPK2dev_dE,
                                                                          d2PK2dev_dEdE,
                                                                          initialCompliance,
                                                                          dPK2dev,
                                                                          timeIncrement.dT,
                                                                          maxwellProperties,
                                                                          stateLayout.getPtr( response.stateVars,
                                                                                              "creepStateVars" ) );
    // compute total PK2 stress
    const Tensor33d PK2 = evaluate( PK2vol + PK2dev );

    // compute Kirchhoff stress from PK2 stress
    const auto [tau,
                dTau_dPK2,
                dTau_dF] = Marmot::ContinuumMechanics::StressMeasures::FirstOrder::KirchhoffStressFromPK2( PK2, F );
    response.tau         = tau;
    response.elasticEnergyDensity = psi_; // TODO: replace by actual elastic energy density

    // compute tangent operator dtau / dF
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, evaluate( dPK2vol_dE + dPK2dev_dE ) ),
                                             0.5 * dC_dF ) +
                       dTau_dF;
  }
} // namespace Marmot::Materials
