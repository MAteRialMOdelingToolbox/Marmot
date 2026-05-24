#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEigenSystems.h"
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

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  double FiniteStrainIsotropicBiotViscoelasticity::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < ( 2 + 1 + maxwellProperties.nMaxwell * 2 ) ) {
      throw std::runtime_error( "Not enough material properties provided to compute density." );
    }
    return materialProperties[nMaterialProperties - 1];
  }
  FiniteStrainIsotropicBiotViscoelasticity::FiniteStrainIsotropicBiotViscoelasticity( const double* materialProperties,
                                                                                      int           nMaterialProperties,
                                                                                      int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( materialProperties[2],
                                                                                    &materialProperties[3] ) ),

      initialCompliance( makeDual( invertMinorSymmetricFourthOrderTensor( std::get< 2 >(
        ContinuumMechanics::EnergyDensityFunctions::SecondOrderDerived::BiotNeoHooke( evaluate( Spatial3D::I ),
                                                                                      K,
                                                                                      G ) ) ) ) )

  {
    initializeStateLayout();
  }

  void FiniteStrainIsotropicBiotViscoelasticity::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                                  const DeformationAD< 3 >&    deformation,
                                                                  const TimeIncrement&         timeIncrement ) const
  {

    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;

    // compute Cauchy-Green deformation
    const Tensor33t< scalar > C = DeformationMeasures::rightCauchyGreen( F );

    // compute eigenvalues and eigenvectors of C
    auto [lam, Q] = Math::computeEigenSystemJacobi( C );
    Tensor33t< scalar > principalStretch( 0. );
    for ( int i = 0; i < 3; ++i ) {
      principalStretch( i, i ) = sqrt( lam( i ) );
    }

    // compute the right stretch tensor U
    const Tensor33t< scalar > U = Q % principalStretch % transpose( Q );

    // compute Biot stress and its derivative with respect to the right stretch tensor
    auto [psi, S_biot, dS_biot_dU] = EnergyDensityFunctions::SecondOrderDerived::BiotNeoHooke( U, K, G );

    // retrieve Biot stress from the previous increment
    Tensor33d& S_biot_old = stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );

    // compute the increment in Biot stress
    const Tensor33t< scalar > dS_biot = S_biot - makeDual( S_biot_old );
    memcpy( S_biot_old.data(), makeReal( S_biot ).data(), 9 * sizeof( double ) );

    // add viscoelastic contribution to Biot stress using the generalized Maxwell model
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel<
      scalar >( S_biot,
                dS_biot_dU,
                initialCompliance,
                dS_biot,
                timeIncrement.dT,
                maxwellProperties,
                stateLayout.getPtr( response.stateVars, "creepStateVars" ) );

    // compute the second Piola-Kirchhoff stress from the Biot stress
    Tensor33t< scalar > S_biot_rotated = transpose( Q ) % S_biot % Q;
    Tensor33t< scalar > PK2_rotated( 0 );

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        PK2_rotated( i, j ) = 2.0 * S_biot_rotated( i, j ) / ( principalStretch( i, i ) + principalStretch( j, j ) );
      }
    }
    Tensor33t< scalar > PK2 = Q % PK2_rotated % transpose( Q );

    // push forward the second Piola-Kirchhoff stress to get the Kirchhoff stress
    const Tensor33t< scalar > tau = StressMeasures::KirchhoffStressFromPK2( PK2, F );
    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi );
  }
} // namespace Marmot::Materials
