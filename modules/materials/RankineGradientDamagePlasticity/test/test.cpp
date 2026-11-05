#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/RankineDamagePlasticity.h"
#include "Marmot/RankineGradientDamagePlasticity.h"
#include <algorithm>
#include <memory>
#include <stdexcept>

using namespace Marmot;
using namespace Testing;
using namespace Materials;

namespace {

  const double E        = 30000.;
  const double nu       = 0.2;
  const double ft       = 3.;
  const double epsilonF = 5e-4;
  const double l        = 2.;
  const double G        = E / ( 2. * ( 1. + nu ) );
  const double lambda   = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );

  // E, nu, ft, epsilonF, nonlocal length, m
  std::vector< double > propertiesWith( double m )
  {
    return { E, nu, ft, epsilonF, l, m, 7.85e-9, 1e-3 };
  }
  std::vector< double > localProperties = { E, nu, ft, epsilonF };

  using Material = RankineGradientDamagePlasticity;

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
    Math::Vector6d                stress    = Math::Vector6d::Zero();
    Eigen::Matrix< double, 5, 1 > stateVars = Eigen::Matrix< double, 5, 1 >::Zero();

    double kappaD() const { return stateVars( 0 ); }
    double acYield() const { return stateVars( 1 ); }
    double omega() const { return stateVars( 2 ); }
    double kappaDWeightedMax() const { return stateVars( 4 ); }
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

  struct LocalResult {
    Math::Vector6d  stress;
    Eigen::Vector3d stateVars;
    Math::Matrix6d  tangent;
    double          dPsi, dDissipation;
  };

  LocalResult updateLocal( const RankineDamagePlasticity& material,
                           const Math::Vector6d&          stressOld,
                           const Eigen::Vector3d&         stateVarsOld,
                           const Math::Vector6d&          dE )
  {
    LocalResult                  r{ stressOld, stateVarsOld, Math::Matrix6d::Zero(), 0., 0. };
    MaterialHypoElastic::state3D state( stressOld, 0., 0., r.stateVars.data() );
    material.computeStress( state, r.tangent, dE, { 0.0, 1.0 } );
    r.stress       = state.stress;
    r.dPsi         = state.elasticEnergyDensity;
    r.dDissipation = state.dissipation;
    return r;
  }

  // compare all tangent blocks [dS/dE dS/dK; dKLocal/dE dKLocal/dK] against central differences
  void checkTangent( const Material& material, const MaterialPointState& stateOld, const Math::Vector6d& dE, double K )
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
    const double          relError = std::max( { error.topLeftCorner( 6, 6 ).cwiseAbs().maxCoeff() / E,
                                                 error.topRightCorner( 6, 1 ).cwiseAbs().maxCoeff() /
                                                   std::max( tangent.topRightCorner( 6, 1 ).norm(), ft / epsilonF ),
                                                 error.bottomRows( 1 ).cwiseAbs().maxCoeff() } );
    throwExceptionOnFailure( relError < 1e-6,
                             Utilities::MakeString() << "tangent check failed, relative error " << relError );
  }

  const std::vector< Math::Vector6d > increments =
    { ( Math::Vector6d() << 1e-4, 0, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 6e-5, 6e-5, 0, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 5e-5, 5e-5, 5e-5, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 6e-5, 5e-5, -1e-5, 0, 0, 0 ).finished(),
      ( Math::Vector6d() << 4e-5, 2e-5, -1e-5, 3e-5, -2e-5, 1e-5 ).finished() };

} // namespace

/// With the nonlocal variable equal to the local one (the limit of a vanishing nonlocal length, where the Helmholtz
/// equation gives K = kappaD), the model reproduces RankineDamagePlasticity for any m, including the total tangent.
void testLocalLimit()
{
  const RankineDamagePlasticity local( localProperties.data(), static_cast< int >( localProperties.size() ), 1 );

  for ( const double m : { 0., 0.5, 1., 2. } ) {
    const auto     props = propertiesWith( m );
    const Material material( props.data(), static_cast< int >( props.size() ), 1 );

    for ( const auto& dE : increments ) {
      MaterialPointState state;
      Math::Vector6d     stressLocal    = Math::Vector6d::Zero();
      Eigen::Vector3d    stateVarsLocal = Eigen::Vector3d::Zero();

      for ( int i = 0; i < 6; i++ ) {
        const auto rl  = updateLocal( local, stressLocal, stateVarsLocal, dE );
        stressLocal    = rl.stress;
        stateVarsLocal = rl.stateVars;

        const auto r = update( material, state, dE, rl.stateVars( 0 ) );
        state        = r.state;

        const Math::Matrix6d totalTangent = r.tan.dStressddStrain + r.tan.dStressddK * r.tan.dKLocalddStrain;
        throwExceptionOnFailure( checkIfEqual< double >( state.stress, rl.stress, 1e-12 ),
                                 Utilities::MakeString() << "local limit (m = " << m << "): stress differs" );
        throwExceptionOnFailure( checkIfEqual< double >( totalTangent, rl.tangent, 1e-9 ),
                                 Utilities::MakeString() << "local limit (m = " << m << "): tangent differs" );
        throwExceptionOnFailure( checkIfEqual( state.kappaD(), rl.stateVars( 0 ), 1e-15 ) &&
                                   checkIfEqual( r.KLocal, rl.stateVars( 0 ), 1e-15 ) &&
                                   checkIfEqual( state.omega(), rl.stateVars( 2 ), 1e-14 ) &&
                                   state.acYield() == rl.stateVars( 1 ),
                                 Utilities::MakeString() << "local limit (m = " << m << "): state differs" );
        throwExceptionOnFailure( checkIfEqual( r.dPsi, rl.dPsi, 1e-10 ) &&
                                   checkIfEqual( r.dDissipation, rl.dDissipation, 1e-10 ),
                                 Utilities::MakeString() << "local limit (m = " << m << "): energies differ" );
        throwExceptionOnFailure( checkIfEqual( r.c, l * l, 1e-15 ), "local limit: c != l^2" );
      }
      throwExceptionOnFailure( state.omega() > 0.1, "local limit: no damage occurred" );
    }
  }

  // m = 0: purely local damage for any nonlocal field, no coupling
  const auto     props = propertiesWith( 0. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );
  const auto     rl = updateLocal( local, Math::Vector6d::Zero(), Eigen::Vector3d::Zero(), increments[0] * 5 );
  const auto     r  = update( material, MaterialPointState(), increments[0] * 5, 123. );
  throwExceptionOnFailure( checkIfEqual< double >( r.state.stress, rl.stress, 1e-12 ) && r.tan.dStressddK.isZero(),
                           "m = 0: not the local model" );
}

/// Uniaxial strain from the virgin state with a given nonlocal variable: closed form solution.
void testUniaxialStrain()
{
  for ( const double m : { 1., 2. } ) {
    const auto     props = propertiesWith( m );
    const Material material( props.data(), static_cast< int >( props.size() ), 1 );

    const double eps = 4e-4;
    const double K   = 2e-4;
    const auto   r = update( material, MaterialPointState(), ( Math::Vector6d() << eps, 0, 0, 0, 0, 0 ).finished(), K );
    const double dLambda = eps - ft / ( lambda + 2. * G );
    const double omega   = 1. - std::exp( -( m * K + ( 1. - m ) * dLambda ) / epsilonF );
    const Math::Vector6d
      target = ( 1. - omega ) *
               ( Math::Vector6d() << ft, lambda * ( eps - dLambda ), lambda * ( eps - dLambda ), 0, 0, 0 ).finished();
    throwExceptionOnFailure( checkIfEqual( r.KLocal, dLambda, 1e-15 ) && checkIfEqual( r.state.omega(), omega, 1e-14 ),
                             "uniaxial strain: damage mismatch" );
    throwExceptionOnFailure( checkIfEqual< double >( r.state.stress, target, 1e-11 ),
                             "uniaxial strain: comparison with analytical solution failed" );
    throwExceptionOnFailure( checkIfEqual< double >( r.tan.dStressddK, -m / epsilonF * target, 1e-9 ),
                             "uniaxial strain: dStress/dK mismatch" );
  }
}

/// Consistent tangents incl. the coupling blocks: one, two and three active surfaces and general loading, with the
/// damage growing, for the standard (m = 1) and the over-nonlocal (m = 2) formulation.
void testTangents()
{
  for ( const double m : { 1., 2. } ) {
    const auto     props = propertiesWith( m );
    const Material material( props.data(), static_cast< int >( props.size() ), 1 );

    for ( const auto& dE : increments ) {
      MaterialPointState state;
      double             K = 0;
      for ( int i = 0; i < 5; i++ ) {
        // nonlocal variable trailing the local one by a bit, so that the damage grows in each increment
        const double KNew = 0.9 * update( material, state, dE, 0. ).KLocal + 1e-5;
        checkTangent( material, state, dE, KNew );
        const auto r = update( material, state, dE, KNew );
        throwExceptionOnFailure( r.state.kappaDWeightedMax() > state.kappaDWeightedMax() || r.KLocal == 0.,
                                 "tangents: damage did not grow" );
        state = r.state;
        K     = KNew;
      }
      throwExceptionOnFailure( state.omega() > 0.05 && state.acYield() == 1., "tangents: no damage occurred" );

      // damage growth by the nonlocal field alone during elastic unloading of the effective stress
      checkTangent( material, state, -0.05 * dE, K + 2e-4 );
      const auto r = update( material, state, -0.05 * dE, K + 2e-4 );
      throwExceptionOnFailure( r.state.acYield() == 0. && r.state.omega() > state.omega(),
                               "tangents: no damage growth by the nonlocal field" );
      throwExceptionOnFailure( r.tan.dKLocalddStrain.isZero(), "tangents: dKLocal/dE nonzero in an elastic step" );

      // no damage growth: decreasing nonlocal field, damaged elastic tangent without coupling
      checkTangent( material, state, -0.05 * dE, 0.5 * K );
      const auto u = update( material, state, -0.05 * dE, 0.5 * K );
      throwExceptionOnFailure( u.state.omega() == state.omega() && u.tan.dStressddK.isZero(),
                               "tangents: damage changed while unloading" );
      const Math::Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      throwExceptionOnFailure( checkIfEqual< double >( u.tan.dStressddStrain, ( 1. - state.omega() ) * Cel, 1e-12 ),
                               "tangents: unloading tangent is not the damaged elastic stiffness" );
    }

    // active set correction: both trial stresses exceed ft, only surface 1 is active at the solution
    checkTangent( material,
                  MaterialPointState(),
                  ( Math::Vector6d() << 1.036e-4, 0.656e-4, 0, 0, 0, 0 ).finished(),
                  3e-5 );
  }
}

/// Zero strain increment: damage may still grow by the nonlocal field.
void testZeroStrainIncrement()
{
  const auto     props = propertiesWith( 1. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );
  // plastic loading, then elastic unloading of the effective stress, so that the tangent is smooth at dE = 0
  auto state = update( material, MaterialPointState(), increments[0] * 3, 1e-4 ).state;
  state      = update( material, state, -0.1 * increments[0], 1e-4 ).state;
  checkTangent( material, state, Math::Vector6d::Zero(), 3e-4 );
  const auto r = update( material, state, Math::Vector6d::Zero(), 3e-4 );
  throwExceptionOnFailure( checkIfEqual( r.state.omega(), 1. - std::exp( -3e-4 / epsilonF ), 1e-14 ),
                           "zero increment: damage not driven by the nonlocal field" );
  throwExceptionOnFailure( checkIfEqual< double >( r.state.stress,
                                                   ( 1. - r.state.omega() ) / ( 1. - state.omega() ) * state.stress,
                                                   1e-12 ),
                           "zero increment: effective stress changed" );
}

void testPropertiesAndFactory()
{
  const auto                                                  props = propertiesWith( 2. );
  std::unique_ptr< MaterialGradientEnhancedHypoElastic< 1 > > material(
    Registry::MaterialGradientEnhancedHypoElasticFactory< 1 >::createMaterial( "RANKINEGRADIENTDAMAGEPLASTICITY",
                                                                               props.data(),
                                                                               static_cast< int >( props.size() ),
                                                                               1 ) );
  throwExceptionOnFailure( material->getNumberOfRequiredStateVars() == 5, "factory: wrong number of state variables" );
  throwExceptionOnFailure( material->getDensity( nullptr ) == props[6], "factory: wrong density" );
  throwExceptionOnFailure( material->getNonlocalViscosity( nullptr )[0] == props[7], "factory: wrong viscosity" );

  const Material withoutOptional( props.data(), 6, 1 );
  throwExceptionOnFailure( withoutOptional.getNonlocalViscosity( nullptr )[0] == 0.,
                           "factory: nonlocal viscosity default is not zero" );

  // invalid properties
  for ( const auto& p : { std::vector< double >{ E, nu, ft, epsilonF, l } } )
    throwExceptionOnFailure( throws< std::runtime_error >(
                               [&] { Material( p.data(), static_cast< int >( p.size() ), 1 ); } ),
                             "properties: missing properties not detected" );
  for ( const auto& p :
        { std::vector< double >{ E, nu, ft, 0., l, 1. }, std::vector< double >{ E, nu, ft, epsilonF, -1., 1. } } )
    throwExceptionOnFailure( throws< std::invalid_argument >(
                               [&] { Material( p.data(), static_cast< int >( p.size() ), 1 ); } ),
                             "properties: nonpositive epsilonF or negative nonlocal length not detected" );
  throwExceptionOnFailure( throws< std::runtime_error >( [&] { withoutOptional.getDensity( nullptr ); } ),
                           "properties: missing density not detected" );
}

/// Element stiffness of GeneralGradientEnhancedDisplacementFiniteElement (C3D8 and C3D20R with the nonlocal field)
/// with this material against central differences of the residual, from a damaged state.
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
  element.assignProperty( Discretizations::MaterialSection( "RANKINEGRADIENTDAMAGEPLASTICITY",
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
  double omegaMax = 0.;
  for ( int q = 0; q < element.getNumberOfQuadraturePoints(); q++ )
    omegaMax = std::max( omegaMax, element.getStateView( "omega", q ).stateLocation[0] );
  throwExceptionOnFailure( omegaMax > 0.1, "element stiffness: no damage" );
}

void testElementStiffness()
{
  for ( const double m : { 1., 2. } ) {
    checkElementStiffness< 8 >( m );
    checkElementStiffness< 20 >( m );
  }
}

/// Complete damage (only reached by a diverging iterate) requests a cutback instead of returning zero stress.
void testCompleteDamageThrows()
{
  const auto     props = propertiesWith( 1. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );
  bool           threw = false;
  try {
    update( material, MaterialPointState(), increments[0], 1e6 );
  }
  catch ( const Utilities::StressUpdateFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "complete damage did not throw" );

  // a completely damaged old state has no effective stress
  MaterialPointState completelyDamaged;
  completelyDamaged.stateVars << 1e6, 1., 1., 1e6, 1e6;
  throwExceptionOnFailure( throws< Utilities::StressUpdateFailed >(
                             [&] { update( material, completelyDamaged, Math::Vector6d::Zero(), 1e6 ); } ),
                           "complete damage of the old state did not throw" );
}

/// Strain increments of a diverged iterate (about 1e288, from a dissipation-controlled run at the onset of
/// dissipation) overflow the energies (inf, and 0 * inf = NaN in their tangents): the material must request a cutback
/// instead of returning non-finite results; a regular increment remains valid.
void testDivergedIterateRequestsCutback()
{
  const auto     props = propertiesWith( 1. );
  const Material material( props.data(), static_cast< int >( props.size() ), 1 );

  const Math::Vector6d dStrainDiverged = ( Math::Vector6d() << -2.1304181664190551e+287,
                                           -1.0379772584363431e+287,
                                           1.1922728518060270e+045,
                                           -2.9804698180504283e+287,
                                           -3.2636090085991650e+287,
                                           -2.4320172494955940e+287 )
                                           .finished();
  const MaterialPointState stateOld;
  bool                     cutback = false;
  try {
    update( material, stateOld, dStrainDiverged, 0. );
  }
  catch ( const Utilities::StressUpdateFailed& ) {
    cutback = true;
  }
  throwExceptionOnFailure( cutback, "diverged iterate did not request a cutback" );

  const auto r = update( material, stateOld, Math::Vector6d::Constant( 1e-5 ), 0. );
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
                                                         testCompleteDamageThrows,
                                                         testDivergedIterateRequestsCutback };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
