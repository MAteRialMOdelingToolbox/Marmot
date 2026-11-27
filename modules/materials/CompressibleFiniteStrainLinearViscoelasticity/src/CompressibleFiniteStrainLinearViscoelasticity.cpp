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
#include <stdexcept>
#include <tuple>
#include <utility>

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

  namespace {

    /* Minor-symmetric part of a fourth-order tensor, 1/4 ( H_ijkl + H_jikl + H_ijlk + H_jilk ). */
    template < typename T >
    Tensor3333t< T > minorSymmetricPart( const Tensor3333t< T >& H )
    {
      Tensor3333t< T > HSym;
      for ( int i = 0; i < 3; i++ )
        for ( int j = i; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int l = k; l < 3; l++ ) {
              const T symmetrized = .25 * ( H( i, j, k, l ) + H( j, i, k, l ) + H( i, j, l, k ) + H( j, i, l, k ) );
              HSym( i, j, k, l )  = symmetrized;
              HSym( j, i, k, l )  = symmetrized;
              HSym( i, j, l, k )  = symmetrized;
              HSym( j, i, l, k )  = symmetrized;
            }
      return HSym;
    }

    /* Energy density and its first three derivatives w.r.t. a symmetric C from a closed-form SecondOrder potential
     * (energy, dPsi/dC, d2Psi/dCdC; generic in the scalar type). The third derivative is the derivative of the
     * closed-form second derivative by first-order dual numbers, which takes 6 evaluations (one per independent
     * entry of C) instead of the 6 x 81 evaluations of third-order dual numbers on the plain potential.
     *
     * The potentials are functions of a general (not necessarily symmetric) C, so their derivatives w.r.t. the
     * independent entries C_ij and C_ji differ (e.g. tr(C C) gives d2/dC01dC01 = 0 but d2/dC01dC10 = 2). Only the
     * symmetric part acts on symmetric increments of C: all derivatives are those of Psi( sym( C ) ), i.e. projected
     * onto their minor-symmetric parts. */
    template < typename SecondOrderPotential >
    std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > differentiateSecondOrderPotential(
      const SecondOrderPotential& potential,
      const Tensor33d&            C )
    {
      // the value parts of every dual evaluation are the energy density and its first derivative at C
      double    psi = 0.;
      Tensor33d dPsi_dC;

      const auto [d2Psi_dCdC, d3Psi_dCdCdC] = Marmot::Math::Differentiation::Automatic::dF_dSymmetricT(
        [&potential, &psi, &dPsi_dC]( const Tensor33t< autodiff::dual >& C_ad ) {
          const auto [psi_ad, dPsi_dC_ad, d2Psi_dCdC_ad] = potential( C_ad );
          psi                                            = autodiff::val( psi_ad );
          for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
              dPsi_dC( i, j ) = .5 * ( autodiff::val( dPsi_dC_ad( i, j ) ) + autodiff::val( dPsi_dC_ad( j, i ) ) );
          return minorSymmetricPart( d2Psi_dCdC_ad );
        },
        C );

      return { psi, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC };
    }
  } // namespace

  std::tuple< double,
              Marmot::Math::Tensors::FastorTensors::Tensor33d,
              Marmot::Math::Tensors::FastorTensors::Tensor3333d,
              Marmot::Math::Tensors::FastorTensors::Tensor333333d >
  CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives(
    const Marmot::Math::Tensors::FastorTensors::Tensor33d& C ) const
  {

    // every base uses its full compressible potential with the parameters in the order of its interface
    // -- see the class docs.
    namespace EDF   = ContinuumMechanics::Elasticity::EnergyDensityFunctions;
    const auto& par = elasticProperties;

    switch ( hyperelasticBase ) {
    case NeoHooke:
      return differentiateSecondOrderPotential(
        [K = par( 0 ), G = par( 1 )]( const auto& C_ad ) {
          return EDF::SecondOrder::CoupledPotential( C_ad, EDF::Coupled::StandardNeoHooke{ K, G } );
        },
        C );
    case Yeoh:
      return differentiateSecondOrderPotential(
        [C1 = par( 0 ), C2 = par( 1 ), C3 = par( 2 ), D1 = par( 3 )]( const auto& C_ad ) {
          return EDF::SecondOrder::YeohPotential( C_ad, C1, C2, C3, D1 );
        },
        C );
    case PenceGouNeoHooke:
      return differentiateSecondOrderPotential( [K = par( 0 ), G = par( 1 )](
                                                  const auto&
                                                    C_ad ) { return EDF::SecondOrder::PenceGouPotentialB( C_ad, K, G ); },
                                                C );
    case MooneyRivlin:
      return differentiateSecondOrderPotential(
        [C1 = par( 0 ), C2 = par( 1 ), D1 = par( 2 )]( const auto& C_ad ) {
          return EDF::SecondOrder::MooneyRivlinPotential( C_ad, C1, C2, D1 );
        },
        C );
    case ArrudaBoyce:
      return differentiateSecondOrderPotential(
        [mu = par( 0 ), lambdaL = par( 1 ), D1 = par( 2 )]( const auto& C_ad ) {
          return EDF::SecondOrder::ArrudaBoyce8ChainPotential( C_ad, mu, lambdaL, D1 );
        },
        C );
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives: "
                                "Unknown hyperelastic base." );
    }
  }

  double CompressibleFiniteStrainLinearViscoelasticity::initialShearModulus() const
  {
    // G0 = 2 ( dPsi/dIbar1 + dPsi/dIbar2 ) at C = I
    const auto& par = elasticProperties;

    switch ( hyperelasticBase ) {
    case NeoHooke:
    case PenceGouNeoHooke: return par( 1 );
    case Yeoh: return 2. * par( 0 );
    case MooneyRivlin: return 2. * ( par( 0 ) + par( 1 ) );
    case ArrudaBoyce: {
      // tends to mu for lambdaL -> infinity
      const double lambdaL2 = par( 1 ) * par( 1 );
      return par( 0 ) * ( 3. * lambdaL2 - 1. ) / ( 3. * ( lambdaL2 - 1. ) );
    }
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::initialShearModulus: "
                                "Unknown hyperelastic base." );
    }
  }

  double CompressibleFiniteStrainLinearViscoelasticity::initialBulkModulus() const
  {
    // K0 = d2U/dJ2 at J = 1
    const auto& par = elasticProperties;

    switch ( hyperelasticBase ) {
    case NeoHooke:
    case PenceGouNeoHooke: return par( 0 );
    case Yeoh: return 2. / par( 3 );
    case MooneyRivlin:
    case ArrudaBoyce: return 2. / par( 2 );
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::initialBulkModulus: "
                                "Unknown hyperelastic base." );
    }
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

    // the initial stiffness dPK2/dE at C = I is isotropic: K0 I(x)I + 2 G0 Dev, hence the compliance is
    // 1/(9 K0) I(x)I + 1/(2 G0) Dev; with shear creep only, the deviatoric part alone is needed
    if ( onlyShearCreep == 1.0 )
      initialCompliance = .5 / initialShearModulus() * Spatial3D::I4;
    else
      initialCompliance = .5 / initialShearModulus() * Spatial3D::DeviatoricSymmetric +
                          1. / ( 9. * initialBulkModulus() ) * Spatial3D::IHyd;

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
    const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

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

    // compute Kirchhoff stress from PK2 stress and tangent operator dtau / dF, with dPK2/dC = 1/2 dPK2/dE
    const Tensor3333d dPK2_dC                  = 0.5 * ( dPK2vol_dE + dPK2dev_dE );
    std::tie( response.tau, tangents.dTau_dF ) = Marmot::ContinuumMechanics::StressMeasures::FirstOrder::
      KirchhoffStressAndTangentFromPK2( PK2, dPK2_dC, F );
    response.elasticEnergyDensity = psi_; // TODO: replace by actual elastic energy density
  }
} // namespace Marmot::Materials
