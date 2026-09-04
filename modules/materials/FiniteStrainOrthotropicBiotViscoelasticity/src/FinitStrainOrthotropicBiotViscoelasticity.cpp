#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace TensorUtility::FastorTensors;
  using namespace TensorUtility::FastorTensors::Indices;
  using namespace TensorUtility::FastorTensors::StandardTensors;

  double FiniteStrainOrthotropicBiotViscoelasticity::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < ( 9 + 1 + maxwellProperties.nMaxwell * 2 ) ) {
      throw std::runtime_error( "Not enough material properties provided to compute density." );
    }
    return materialProperties[nMaterialProperties - 1];
  }

  FiniteStrainOrthotropicBiotViscoelasticity::FiniteStrainOrthotropicBiotViscoelasticity(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel ),
      E1( materialProperties[0] ),
      E2( materialProperties[1] ),
      E3( materialProperties[2] ),
      nu12( materialProperties[3] ),
      nu13( materialProperties[4] ),
      nu23( materialProperties[5] ),
      G12( materialProperties[6] ),
      G13( materialProperties[7] ),
      G23( materialProperties[8] ),
      maxwellProperties(
        ContinuumMechanics::Viscoelasticity::FiniteStrain::createMaxwellProperties( materialProperties[9],
                                                                                    &materialProperties[10] ) ),
      dBiotStress_dU( makeDual( ContinuumMechanics::Voigt::voigtToStiffnessFastor(
        ContinuumMechanics::Elasticity::Orthotropic::
          stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 ) ) ) )
  {

    initializeStateLayout();
  }

  void FiniteStrainOrthotropicBiotViscoelasticity::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                                    const DeformationAD< 3 >&    deformation,
                                                                    const TimeIncrement&         timeIncrement ) const
  {

    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;
    using namespace ContinuumMechanics::Kinematics;
    // compute Cauchy-Green deformation
    const Tensor33t< scalar > C = DeformationMeasures::rightCauchyGreen( F );

    // compute principal stretches and directions
    auto [lam, Q] = TensorUtility::FastorTensors::StandardTensors::computeEigenSystemJacobi( C );
    Tensor33t< scalar > principalStretch( 0. );
    for ( int i = 0; i < 3; ++i ) {
      principalStretch( i, i ) = sqrt( lam( i ) );
    }

    // compute right stretch tensor in the original configuration
    const Tensor33t< scalar > U = Q % principalStretch % transpose( Q );

    const Tensor33t< scalar > I = makeDual( Spatial3D::I );

    // compute Biot stress
    Tensor33t< scalar > S_biot = einsum< ijkl, kl >( dBiotStress_dU, evaluate( U - I ) );

    // get old Biot stress from state variables
    Tensor33d S_biot_old( stateLayout.getPtr( response.stateVars, "S0_old" ) );

    // compute change in Biot stress
    const Tensor33t< scalar > dS_biot = S_biot - makeDual( S_biot_old );

    // store current Biot stress in state variables
    memcpy( stateLayout.getPtr( response.stateVars, "S0_old" ), makeReal( S_biot ).data(), 9 * sizeof( double ) );

    // add viscoelastic contribution to Biot stress
    ContinuumMechanics::Viscoelasticity::FiniteStrain::evaluateGeneralizedMaxwellModel<
      scalar >( S_biot,
                dS_biot,
                timeIncrement.dT,
                maxwellProperties,
                stateLayout.getPtr( response.stateVars, "creepStateVars" ) );

    // transform to PK2 stress
    Tensor33t< scalar > S_biot_rotated = transpose( Q ) % S_biot % Q;
    Tensor33t< scalar > PK2_rotated( 0 );

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        PK2_rotated( i, j ) = 2.0 * S_biot_rotated( i, j ) / ( principalStretch( i, i ) + principalStretch( j, j ) );
      }
    }
    Tensor33t< scalar > PK2 = Q % PK2_rotated % transpose( Q );

    // transform to Kirchhoff stress
    const Tensor33t< scalar > tau = StressMeasures::KirchhoffStressFromPK2( PK2, F );
    response.tau                  = tau;
    response.elasticEnergyDensity = 0.5 *
                                    einsum< ij, ij >( makeReal( S_biot ), makeReal( evaluate( U - I ) ) ).toscalar();
  }
} // namespace Marmot::Materials
