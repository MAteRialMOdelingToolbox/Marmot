#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/RankineGradientSofteningPlasticity.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>

using namespace Marmot;
using namespace Testing;
using namespace Materials;

namespace {

  const double E              = 30000.;
  const double nu             = 0.2;
  const double fy0            = 3.;
  const double HLin           = -2000.;
  const double nonlocalLength = 2.;
  const double G              = E / ( 2. * ( 1. + nu ) );
  const double lambda         = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

  // E, nu, fy0, HLin, deltaYieldStress, delta, nonlocal length, m, density, nonlocal viscosity
  std::vector< double > propertiesWith( double m, bool exponential = false )
  {
    return { E,
             nu,
             fy0,
             exponential ? -500. : HLin,
             exponential ? -2. : 0.,
             exponential ? 3000. : 0.,
             nonlocalLength,
             m,
             7.85e-9,
             1e-3 };
  }

  using Material = RankineGradientSofteningPlasticity;

  // true if f throws an exception of type Ex
  template < typename Ex, typename F >
  bool throws( F&& f )
  {
    try {
      f();
    }
    catch ( const Ex& ) {
      return true;
    }
    return false;
  }

  struct MaterialPointState {
    Math::Vector6d  stress    = Math::Vector6d::Zero();
    Eigen::Vector4d stateVars = Eigen::Vector4d::Zero(); // kappa, acYield, kappaNonlocal, kappaWeighted

    double kappa() const { return stateVars( 0 ); }
    double acYield() const { return stateVars( 1 ); }
    double kappaWeighted() const { return stateVars( 3 ); }
  };

  struct UpdateResult {
    MaterialPointState state;
    Material::tangents tan;
    double             KLocal, c, dPsi, dDissipation;
  };

  UpdateResult update( const Material&           material,
                       const MaterialPointState& stateOld,
                       const Math::Vector6d&     dE,
                       double                    K )
  {
    UpdateResult       result{ .state = stateOld, .tan = {}, .KLocal = 0., .c = 0., .dPsi = 0., .dDissipation = 0. };
    Material::response res{ .stress               = stateOld.stress,
                            .KLocal               = Eigen::Vector< double, 1 >::Zero(),
                            .c                    = Eigen::Vector< double, 1 >::Zero(),
                            .stateVars            = result.state.stateVars.data(),
                            .elasticEnergyDensity = 0.,
                            .dissipation          = 0. };
    const Material::increment inc{ .dStrain = dE,
                                   .K       = Eigen::Vector< double, 1 >::Constant( K ),
                                   .dK      = Eigen::Vector< double, 1 >::Zero(),
                                   .time    = 0.,
                                   .dT      = 1. };
    material.computeStress( res, result.tan, inc );
    result.state.stress = res.stress;
    result.KLocal       = res.KLocal( 0 );
    result.c            = res.c( 0 );
    result.dPsi         = res.elasticEnergyDensity;
    result.dDissipation = res.dissipation;
    return result;
  }

  // compare all tangent blocks [dS/dE dS/dK; dKLocal/dE dKLocal/dK] against central differences
  void checkTangent( const Material&           material,
                     const MaterialPointState& stateOld,
                     const Math::Vector6d&     dE,
                     double                    K,
                     const std::string&        label )
  {
    using namespace Math::Differentiation::Numerical;
    const auto r = update( material, stateOld, dE, K );

    Eigen::Matrix< double, 7, 7 > tangent;
    tangent << r.tan.dStressddStrain, r.tan.dStressddK, r.tan.dKLocalddStrain, r.tan.dKLocalddK;

    // scaled arguments, since the default step is relative to max(1, |x|)
    const double          scale = 1e-4;
    const Eigen::VectorXd x     = ( Eigen::VectorXd( 7 ) << dE, K ).finished() / scale;
    const Eigen::MatrixXd
      numericalTangent = centralDifference(
                           [&]( const Eigen::VectorXd& y ) -> Eigen::VectorXd {
                             const auto ry = update( material, stateOld, scale * y.head( 6 ), scale * y( 6 ) );
                             return ( Eigen::VectorXd( 7 ) << ry.state.stress, ry.KLocal ).finished();
                           },
                           x ) /
                         scale;

    const Eigen::MatrixXd error    = tangent - numericalTangent;
    const double          relError = std::max(
      { error.topLeftCorner( 6, 6 ).cwiseAbs().maxCoeff() / E,
                 error.topRightCorner( 6, 1 ).cwiseAbs().maxCoeff() / std::max( tangent.topRightCorner( 6, 1 ).norm(), -HLin ),
                 error.bottomRows( 1 ).cwiseAbs().maxCoeff() } );
    throwExceptionOnFailure( relError < 1e-6,
                             Utilities::MakeString()
                               << label << ": tangent check failed, relative error " << relError );
  }

  // one, two and three active surfaces, and general loading
  const std::vector< Math::Vector6d > increments =
    { ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 6e-5, 6e-5, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 5e-5, 5e-5, 5e-5, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 6e-5, 5e-5, -1e-5, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished() };

} // namespace

/// m = 0: the nonlocal variable has no effect and the model is RankineSofteningPlasticity.
void testLocalLimit()
{
  for ( const bool exponential : { false, true } ) {
    const auto                       props = propertiesWith( 0., exponential );
    const Material                   material( props.data(), static_cast< int >( props.size() ), 1 );
    const RankineSofteningPlasticity local( props.data(), 6, 1 );

    for ( const auto& dE : increments ) {
      MaterialPointState state;
      Math::Vector6d     stressLocal    = Math::Vector6d::Zero();
      Eigen::Vector2d    stateVarsLocal = Eigen::Vector2d::Zero();
      for ( int i = 0; i < 6; i++ ) {
        const auto r = update( material, state, dE, 1e-4 * i );
        state        = r.state;

        Math::Matrix6d               tangentLocal;
        MaterialHypoElastic::state3D stateLocal( stressLocal, 0., 0., stateVarsLocal.data() );
        local.computeStress( stateLocal, tangentLocal, dE, { 0.0, 1.0 } );
        stressLocal = stateLocal.stress;

        throwExceptionOnFailure( checkIfEqual< double >( state.stress, stressLocal, 1e-11 ) &&
                                   checkIfEqual( state.kappa(), stateVarsLocal( 0 ), 1e-15 ) &&
                                   state.acYield() == stateVarsLocal( 1 ),
                                 "local limit: state differs from RankineSofteningPlasticity" );
        throwExceptionOnFailure( checkIfEqual< double >( r.tan.dStressddStrain, tangentLocal, 1e-8 ) &&
                                   r.tan.dStressddK.isZero() && r.tan.dKLocalddK.isZero(),
                                 "local limit: tangents differ from RankineSofteningPlasticity" );
        throwExceptionOnFailure( checkIfEqual( r.dPsi, stateLocal.elasticEnergyDensity, 1e-13 ) &&
                                   checkIfEqual( r.dDissipation, stateLocal.dissipation, 1e-13 ),
                                 "local limit: energies differ from RankineSofteningPlasticity" );
      }
      throwExceptionOnFailure( state.kappa() > 0., "local limit: no plastic flow occurred" );
    }
  }
}

/// Uniaxial strain with linear softening: sigma_1 = fy0 + H ( m K + (1 - m) dLambda ).
void testUniaxialStrain()
{
  const double eps = 2e-4, K = 1e-4;
  for ( const double m : { 0.5, 1., 2. } ) {
    const auto     props = propertiesWith( m );
    const Material material( props.data(), static_cast< int >( props.size() ), 1 );

    const auto   r       = update( material, MaterialPointState(), increments[0] * ( eps / 1e-4 ), K );
    const double dLambda = ( ( lambda + 2. * G ) * eps - fy0 - HLin * m * K ) / ( lambda + 2. * G + HLin * ( 1. - m ) );
    const double fy      = fy0 + HLin * ( m * K + ( 1. - m ) * dLambda );

    Math::Vector6d stress = Math::Vector6d::Zero();
    stress.head< 3 >() << ( lambda + 2. * G ) * ( eps - dLambda ), lambda * ( eps - dLambda ),
      lambda * ( eps - dLambda );

    throwExceptionOnFailure( checkIfEqual( stress( 0 ), fy, 1e-11 ),
                             "uniaxial strain: inconsistent analytical solution" );
    throwExceptionOnFailure( checkIfEqual( r.state.kappa(), dLambda, 1e-15 ) &&
                               checkIfEqual( r.KLocal, dLambda, 1e-15 ) &&
                               checkIfEqual( r.state.kappaWeighted(), m * K + ( 1. - m ) * dLambda, 1e-15 ),
                             "uniaxial strain: softening variables mismatch" );
    throwExceptionOnFailure( checkIfEqual< double >( r.state.stress, stress, 1e-11 ),
                             "uniaxial strain: comparison with analytical solution failed" );
    throwExceptionOnFailure( checkIfEqual( r.c, nonlocalLength * nonlocalLength, 1e-15 ),
                             "uniaxial strain: wrong nonlocal parameter" );
    throwExceptionOnFailure( checkIfEqual( r.dDissipation, fy * dLambda, 1e-13 ),
                             "uniaxial strain: wrong dissipation" );
    // m > 1: local hardening at a fixed nonlocal variable
    if ( m > 1. )
      throwExceptionOnFailure( update( material, MaterialPointState(), increments[0] * ( 2. * eps / 1e-4 ), K )
                                   .state.stress( 0 ) > r.state.stress( 0 ),
                               "uniaxial strain: no local hardening for m > 1" );
  }
}

/// Consistent tangents incl. the coupling blocks, during softening and with the yield stress floored at zero.
void testTangents()
{
  for ( const bool exponential : { false, true } )
    for ( const double m : { 1., 2. } ) {
      const auto     props = propertiesWith( m, exponential );
      const Material material( props.data(), static_cast< int >( props.size() ), 1 );

      for ( const auto& dE : increments ) {
        MaterialPointState state;
        double             K = 0;
        for ( int i = 0; i < 5; i++ ) {
          // nonlocal variable following the local one
          const double KNew = 0.8 * update( material, state, dE, K ).KLocal + 1e-5;
          checkTangent( material, state, dE, KNew, "loading" );
          state = update( material, state, dE, KNew ).state;
          K     = KNew;
        }
        throwExceptionOnFailure( state.acYield() == 1. && state.kappa() > 1e-5, "tangents: no plastic flow occurred" );

        // elastic unloading
        checkTangent( material, state, -0.05 * dE, K, "unloading" );
        const auto u = update( material, state, -0.05 * dE, K );
        throwExceptionOnFailure( u.state.acYield() == 0. && u.tan.dStressddK.isZero() &&
                                   u.tan.dKLocalddStrain.isZero() && u.dDissipation == 0.,
                                 "tangents: plastic flow while unloading" );
      }
    }

  // linear softening with the yield stress floored at zero by the nonlocal variable, triaxial extension (all three
  // surfaces active with positive consistency parameters, away from the corners of the active set): no stress, no
  // sensitivity to the nonlocal variable
  {
    const auto     props = propertiesWith( 1. );
    const Material material( props.data(), static_cast< int >( props.size() ), 1 );
    const double   K = 2. * fy0 / -HLin;
    checkTangent( material, MaterialPointState(), increments[2], K, "floored yield stress" );
    const auto r = update( material, MaterialPointState(), increments[2], K );
    throwExceptionOnFailure( r.state.stress.norm() < 1e-10 && r.tan.dStressddK.isZero() && r.dDissipation == 0.,
                             "floored yield stress: stress or sensitivity to the nonlocal variable remain" );
  }
}

/// Zero strain increment: a growing nonlocal variable lowers the yield stress, which returns the stress.
void testZeroStrainIncrement()
{
  const auto     props = propertiesWith( 2. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );
  const auto     state = update( material, MaterialPointState(), increments[0] * 2, 1e-4 ).state;

  const double KNew = 3e-4;
  checkTangent( material, state, Math::Vector6d::Zero(), KNew, "zero increment" );
  const auto r = update( material, state, Math::Vector6d::Zero(), KNew );
  throwExceptionOnFailure( r.state.acYield() == 1. && r.state.kappa() > state.kappa() &&
                             r.state.stress( 0 ) < state.stress( 0 ),
                           "zero increment: no return by the nonlocal variable" );
  throwExceptionOnFailure( checkIfEqual( r.state.stress( 0 ), fy0 + HLin * r.state.kappaWeighted(), 1e-11 ),
                           "zero increment: stress not on the yield surface" );

  // a decreasing nonlocal variable raises the yield stress: elastic
  const auto e = update( material, state, Math::Vector6d::Zero(), 0.5e-4 );
  throwExceptionOnFailure( e.state.acYield() == 0. && checkIfEqual< double >( e.state.stress, state.stress, 1e-13 ),
                           "zero increment: plastic flow for a decreasing nonlocal variable" );

  // a non-finite nonlocal field requests a cutback; NaN would otherwise floor the yield stress at zero
  for ( const double KNonFinite :
        { std::numeric_limits< double >::quiet_NaN(), std::numeric_limits< double >::infinity() } )
    throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >(
                               [&] { update( material, state, Math::Vector6d::Zero(), KNonFinite ); } ),
                             "non-finite nonlocal field did not request a cutback" );
}

void testPropertiesAndFactory()
{
  const auto                                                  props = propertiesWith( 2. );
  std::unique_ptr< MaterialGradientEnhancedHypoElastic< 1 > > material(
    Registry::MaterialGradientEnhancedHypoElasticFactory< 1 >::createMaterial( "RANKINEGRADIENTSOFTENINGPLASTICITY",
                                                                               props.data(),
                                                                               static_cast< int >( props.size() ),
                                                                               1 ) );
  throwExceptionOnFailure( material->getNumberOfRequiredStateVars() == 4, "factory: wrong number of state variables" );
  throwExceptionOnFailure( material->getDensity( nullptr ) == props[8], "factory: wrong density" );
  throwExceptionOnFailure( material->getNonlocalViscosity( nullptr )[0] == props[9], "factory: wrong viscosity" );

  const Material withoutOptional( props.data(), 8, 1 );
  throwExceptionOnFailure( withoutOptional.getNonlocalViscosity( nullptr )[0] == 0.,
                           "factory: nonlocal viscosity default is not zero" );

  // invalid properties
  throwExceptionOnFailure( throws< std::runtime_error >( [&] { Material( props.data(), 7, 1 ); } ),
                           "properties: missing properties not detected" );
  const std::vector< double > negativeLength = { E, nu, fy0, HLin, 0., 0., -1., 1. };
  throwExceptionOnFailure( throws< std::invalid_argument >( [&] {
                             Material( negativeLength.data(), static_cast< int >( negativeLength.size() ), 1 );
                           } ),
                           "properties: negative nonlocal length not detected" );
  throwExceptionOnFailure( throws< std::runtime_error >( [&] { withoutOptional.getDensity( nullptr ); } ),
                           "properties: missing density not detected" );
}

/// Element stiffness of GeneralGradientEnhancedDisplacementFiniteElement (C3D8 and C3D20R with the nonlocal field)
/// with this material against central differences of the residual, from a plastic state.
template < int nNodes >
void checkElementStiffness( double m )
{
  using Element      = Discretizations::Elements::GeneralGradientEnhancedDisplacementFiniteElement< 3, nNodes >;
  const auto intType = nNodes == 20 ? Discretizations::Elements::Quadrature::IntegrationTypes::ReducedIntegration
                                    : Discretizations::Elements::Quadrature::IntegrationTypes::FullIntegration;
  Element    element( 1, intType, Element::SectionType::Solid );

  const std::vector< Eigen::Vector3d > corners = { { 0, 0, 0 },
                                                   { 1.1, 0, 0.05 },
                                                   { 1, 1, 0 },
                                                   { -0.05, 0.9, 0 },
                                                   { 0, 0.05, 1 },
                                                   { 1, 0, 1.1 },
                                                   { 1.05, 1, 1 },
                                                   { 0, 1, 0.95 } };
  std::vector< Eigen::Vector3d >       nodes   = corners;
  if constexpr ( nNodes == 20 ) {
    const int edges[12][2] = { { 0, 1 },
                               { 1, 2 },
                               { 2, 3 },
                               { 3, 0 },
                               { 4, 5 },
                               { 5, 6 },
                               { 6, 7 },
                               { 7, 4 },
                               { 0, 4 },
                               { 1, 5 },
                               { 2, 6 },
                               { 3, 7 } };
    for ( const auto& e : edges )
      nodes.push_back( 0.5 * ( corners[e[0]] + corners[e[1]] ) );
  }
  std::vector< double > coordinates;
  for ( const auto& x : nodes )
    coordinates.insert( coordinates.end(), x.data(), x.data() + 3 );

  const auto props = propertiesWith( m );
  element.assignNodeCoordinates( coordinates.data() );
  element.assignProperty( Discretizations::MaterialSection( "RANKINEGRADIENTSOFTENINGPLASTICITY",
                                                            props.data(),
                                                            static_cast< int >( props.size() ) ) );
  std::vector< double > stateVars( element.getNumberOfRequiredStateVars(), 0. );
  element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
  element.initializeYourself();

  constexpr int   n = Element::sizeLoadVector, nU = Element::sizeDoFU;
  Eigen::VectorXd dQ1 = Eigen::VectorXd::Zero( n ), dQ2 = Eigen::VectorXd::Zero( n );
  for ( int i = 0; i < nU; i++ )
    dQ1[i] = 1e-4 * std::sin( 1.5 * ( i + 1 ) ) + ( i % 3 == 0 ? 2e-4 * coordinates[i] : 0. );
  for ( int A = 0; A < nNodes; A++ )
    dQ1[nU + A] = 2e-4 * ( 1. + 0.2 * std::cos( 0.9 * ( A + 1 ) ) );
  dQ2 = 0.75 * dQ1;

  auto kernels = [&]( const Eigen::VectorXd& Q, const Eigen::VectorXd& dQ ) {
    Eigen::VectorXd P = Eigen::VectorXd::Zero( n );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( n, n );
    element.computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0., 1. );
    return std::make_pair( P, K );
  };
  kernels( dQ1, dQ1 );
  const std::vector< double > stateVarsOld = stateVars;
  auto                        residual     = [&]( const Eigen::VectorXd& dQ ) {
    std::copy( stateVarsOld.begin(), stateVarsOld.end(), stateVars.begin() );
    return kernels( dQ1 + dQ, dQ );
  };

  const auto [P, K] = residual( dQ2 );
  Eigen::MatrixXd KNumerical( n, n );
  for ( int j = 0; j < n; j++ ) {
    const double    h  = 1e-7 * std::max( std::abs( dQ2[j] ), 1e-4 );
    Eigen::VectorXd dp = dQ2, dm = dQ2;
    dp[j] += h;
    dm[j] -= h;
    KNumerical.col( j ) = ( residual( dp ).first - residual( dm ).first ) / ( 2. * h );
  }
  const auto blockError = [&]( int r0, int c0, int nr, int nc ) {
    return ( K.block( r0, c0, nr, nc ) - KNumerical.block( r0, c0, nr, nc ) ).cwiseAbs().maxCoeff() /
           K.block( r0, c0, nr, nc ).cwiseAbs().maxCoeff();
  };
  const double error = std::max( { blockError( 0, 0, nU, nU ),
                                   blockError( 0, nU, nU, nNodes ),
                                   blockError( nU, 0, nNodes, nU ),
                                   blockError( nU, nU, nNodes, nNodes ) } );
  throwExceptionOnFailure( error < 1e-5,
                           Utilities::MakeString() << "element stiffness differs from central differences, m = " << m
                                                   << ", relative error " << error );
  double kappaMax = 0.;
  for ( int q = 0; q < element.getNumberOfQuadraturePoints(); q++ )
    kappaMax = std::max( kappaMax, element.getStateView( "kappa", q ).stateLocation[0] );
  throwExceptionOnFailure( kappaMax > 1e-5, "element stiffness: no plastic flow" );
}

void testElementStiffness()
{
  for ( const double m : { 1., 2. } ) {
    checkElementStiffness< 8 >( m );
    checkElementStiffness< 20 >( m );
  }
}

/// Strain increments of a diverged iterate (about 1e288) overflow the results: the material must request a cutback
/// instead of returning non-finite results; a regular increment remains valid.
void testDivergedIterateRequestsCutback()
{
  const auto     props = propertiesWith( 2. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );

  const Math::Vector6d dStrainDiverged = ( Math::Vector6d() << -2.1304181664190551e+287,
                                           -1.0379772584363431e+287,
                                           1.1922728518060270e+045,
                                           -2.9804698180504283e+287,
                                           -3.2636090085991650e+287,
                                           -2.4320172494955940e+287 )
                                           .finished();
  throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >(
                             [&] { update( material, MaterialPointState(), dStrainDiverged, 0. ); } ),
                           "diverged iterate did not request a cutback" );

  const auto r = update( material, MaterialPointState(), Math::Vector6d::Constant( 1e-5 ), 0. );
  throwExceptionOnFailure( r.state.stress.allFinite(), "regular increment gave non-finite results" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testLocalLimit,
                                                         testUniaxialStrain,
                                                         testTangents,
                                                         testZeroStrainIncrement,
                                                         testPropertiesAndFactory,
                                                         testElementStiffness,
                                                         testDivergedIterateRequestsCutback };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
