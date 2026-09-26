#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

using Mat = GradientEnhancedFiniteStrainDruckerPrager;

namespace {

  constexpr double K = 3500., G = 1500.;

  // K, G, c0, phi, psi, H, epsF, omegaMax, l, m, rho
  std::array< double, 11 > properties( double c0,
                                       double phi,
                                       double psi,
                                       double H    = 100.,
                                       double epsF = 1e10,
                                       double m    = 0.0 )
  {
    return { K, G, c0, phi, psi, H, epsF, 0.99, 2.0, m, 2.4e-9 };
  }

  struct Result {
    Tensor33d                   tau;
    double                      L;
    double                      dissipation;
    std::vector< double >       state;
    Mat::AlgorithmicModuli< 3 > t;
  };

  // one call on a COPY of the given state (the history is not committed)
  Result evaluate( const Mat&                   mat,
                   const Tensor33d&             F,
                   double                       N,
                   const std::vector< double >& stateIn,
                   double                       dissipationIn = 0.0 )
  {
    std::vector< double >          state = stateIn;
    Mat::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0., 0., 0., dissipationIn, state.data() );
    Mat::AlgorithmicModuli< 3 >    t;
    mat.computeStress( response, t, { F, N }, { 0.0, 1.0 } );
    return { response.tau, response.L, response.dissipation, state, t };
  }

  std::vector< double > freshState( Mat& mat )
  {
    std::vector< double > state( mat.getNumberOfRequiredStateVars() );
    mat.initializeYourself( state.data(), state.size() );
    return state;
  }

  double stateValue( const Mat& mat, const std::vector< double >& state, const std::string& name )
  {
    return *mat.getStateView( name, const_cast< double* >( state.data() ) ).stateLocation;
  }

  Eigen::Matrix3d toEigen( const Tensor33d& T )
  {
    Eigen::Matrix3d M;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        M( i, j ) = T( i, j );
    return M;
  }

  Tensor33d stretch( double lx, double ly, double lz )
  {
    Tensor33d F( 0.0 );
    F( 0, 0 ) = lx;
    F( 1, 1 ) = ly;
    F( 2, 2 ) = lz;
    return F;
  }

  // a general deformation gradient with shear, compression and rotation
  Tensor33d testF( double s )
  {
    Tensor33d F = { { 1. - 0.3 * s, 0.9 * s, -0.1 * s },
                    { 0.2 * s, 1. - 0.1 * s, 0.4 * s },
                    { 0.05 * s, -0.2 * s, 1. } };
    return F;
  }

  double yieldFunction( const Tensor33d& tau, double alphaP, double c0, double phi, double H )
  {
    const Eigen::Vector3d S = Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d >( toEigen( tau ) ).eigenvalues();
    const double          p = S.mean();
    const double          q = std::sqrt( 0.5 * ( S.array() - p ).matrix().squaredNorm() );
    const auto [eta, xi]    = Mat::outerConeParameters( phi );
    return q + eta * p - xi * ( c0 + H * alphaP );
  }

  const std::string where = " in " + std::string( __FILE__ );

} // namespace

void testElasticRangeIsCompressibleNeoHooke()
{
  const auto props = properties( 1e6, 30., 10. ); // never yields
  Mat        mat( props.data(), props.size(), 1 );

  CompressibleNeoHooke                            reference( props.data(), 2, 1 );
  CompressibleNeoHooke::ConstitutiveResponse< 3 > r( Tensor33d( 0.0 ), 0., 0., nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    t;

  for ( double s : { 0.05, 0.2, 0.5 } ) {
    reference.computeStress( r, t, { testF( s ) }, { 0.0, 1.0 } );
    const auto res = evaluate( mat, testF( s ), 0.0, freshState( mat ) );
    throwExceptionOnFailure( checkIfEqual( res.tau, r.tau, 1e-8 * ( 1. + norm( r.tau ) ) ),
                             "elastic stress != CompressibleNeoHooke" + where );
    throwExceptionOnFailure( checkIfEqual( res.t.dTau_dF, t.dTau_dF, 1e-3 * ( 1. + norm( t.dTau_dF ) ) ),
                             "elastic tangent != CompressibleNeoHooke" + where );
  }
}

void testConsistencyAndFlowRule()
{
  const double c0 = 5., phi = 30., H = 100.;

  for ( double psi : { 0.0, 10.0, 30.0 } ) {
    const auto props = properties( c0, phi, psi, H );
    Mat        mat( props.data(), props.size(), 1 );
    const auto res = evaluate( mat, testF( 0.1 ), 0.0, freshState( mat ) );

    const double alphaP = stateValue( mat, res.state, "alphaP" );
    throwExceptionOnFailure( alphaP > 0.0, "the test deformation must yield" + where );

    const double f = yieldFunction( res.tau, alphaP, c0, phi, H );
    throwExceptionOnFailure( std::abs( f ) < 1e-8 * c0,
                             MakeString() << "the returned state is not on the yield surface, f = " << f << where );

    // exponential map: ln det Fp is the volumetric plastic log strain = etaBar dLambda = (etaBar / xi) alphaP
    const Tensor33d Fp( mat.getStateView( "Fp", const_cast< double* >( res.state.data() ) ).stateLocation );
    const double    lnJp      = std::log( determinant( Fp ) );
    const double    etaBar    = Mat::outerConeParameters( psi ).first;
    const double    xi        = Mat::outerConeParameters( phi ).second;
    const double    predicted = etaBar / xi * alphaP;
    throwExceptionOnFailure( std::abs( lnJp - predicted ) < 1e-10,
                             MakeString() << "flow rule: ln det Fp = " << lnJp << ", expected " << predicted << where );
  }
}

void testApexReturn()
{
  const double c0 = 5., phi = 30., psi = 20., H = 100.;
  const auto   props = properties( c0, phi, psi, H );
  Mat          mat( props.data(), props.size(), 1 );

  // hydrostatic tension far beyond the apex
  const auto res = evaluate( mat, stretch( 1.01, 1.01, 1.01 ), 0.0, freshState( mat ) );

  const Eigen::Matrix3d tau = toEigen( res.tau );
  const double          p   = tau.trace() / 3.0;
  throwExceptionOnFailure( ( tau - p * Eigen::Matrix3d::Identity() ).norm() < 1e-9 * std::abs( p ),
                           "the apex state must be hydrostatic" + where );

  const auto [eta, xi] = Mat::outerConeParameters( phi );
  const double alphaP  = stateValue( mat, res.state, "alphaP" );
  throwExceptionOnFailure( std::abs( eta * p - xi * ( c0 + H * alphaP ) ) < 1e-8 * c0,
                           "the apex state must satisfy eta p = xi c" + where );
}

void testObjectivity()
{
  const auto props = properties( 5., 30., 10. );
  Mat        mat( props.data(), props.size(), 1 );

  const auto      ref = evaluate( mat, testF( 0.1 ), 0.0, freshState( mat ) );
  const Tensor33d Fp0( mat.getStateView( "Fp", const_cast< double* >( ref.state.data() ) ).stateLocation );

  for ( int deg = 30; deg <= 180; deg += 50 ) {
    const double phi = Marmot::Math::degToRad( deg );
    Tensor33d    Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.;

    const auto rot = evaluate( mat, Tensor33d( einsum< ik, kj, to_ij >( Q, testF( 0.1 ) ) ), 0.0, freshState( mat ) );
    const Tensor33d tauRotated = einsum< iI, IJ, jJ, to_ij >( Q, ref.tau, Q );
    const Tensor33d Fp( mat.getStateView( "Fp", const_cast< double* >( rot.state.data() ) ).stateLocation );

    throwExceptionOnFailure( checkIfEqual( rot.tau, tauRotated, 1e-9 * norm( ref.tau ) ),
                             "tau is not objective" + where );
    throwExceptionOnFailure( checkIfEqual( Fp, Fp0, 1e-10 ), "Fp must not see a superposed rotation" + where );
  }
}

// the algorithmic tangents against central differences of the full update, from the same initial state
void checkTangents( const Mat&                   mat,
                    const Tensor33d&             F,
                    double                       N,
                    const std::vector< double >& state0,
                    const std::string&           branch )
{
  const auto res = evaluate( mat, F, N, state0 );
  const auto tol = [&]( double num, double scale ) { return 1e-5 * ( scale + std::abs( num ) ); };

  const double h = 1e-6;
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Tensor33d Fp = F, Fm = F;
      Fp( k, l ) += h;
      Fm( k, l ) -= h;
      const auto rp = evaluate( mat, Fp, N, state0 );
      const auto rm = evaluate( mat, Fm, N, state0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          const double num = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2 * h );
          throwExceptionOnFailure( std::abs( res.t.dTau_dF( i, j, k, l ) - num ) < tol( num, 1.0 ),
                                   MakeString() << branch << ": dTau_dF(" << i << j << k << l
                                                << ") = " << res.t.dTau_dF( i, j, k, l ) << " vs " << num << where );
        }
      const double numL = ( rp.L - rm.L ) / ( 2 * h );
      throwExceptionOnFailure( std::abs( res.t.dL_dF( k, l ) - numL ) < tol( numL, 1e-3 ),
                               MakeString() << branch << ": dL_dF(" << k << l << ") = " << res.t.dL_dF( k, l ) << " vs "
                                            << numL << where );
    }

  const double hN = 1e-7;
  const auto   rp = evaluate( mat, F, N + hN, state0 );
  const auto   rm = evaluate( mat, F, N - hN, state0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      const double num = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2 * hN );
      throwExceptionOnFailure( std::abs( res.t.dTau_dN( i, j ) - num ) < tol( num, 1.0 ),
                               branch + ": dTau_dN inconsistent" + where );
    }
  const double numL = ( rp.L - rm.L ) / ( 2 * hN );
  throwExceptionOnFailure( std::abs( res.t.dL_dN - numL ) < 1e-6 * ( 1. + std::abs( numL ) ),
                           branch + ": dL_dN inconsistent" + where );
}

void testTangentInThePlasticBranch()
{
  // plasticity AND damage active
  const auto props = properties( 5., 30., 10., 100., 0.05, 0.5 );
  Mat        mat( props.data(), props.size(), 1 );
  checkTangents( mat, testF( 0.1 ), 2e-3, freshState( mat ), "cone" );

  // from a plastic, damaged state (Fp != I)
  std::vector< double > state = evaluate( mat, testF( 0.1 ), 2e-3, freshState( mat ) ).state;
  checkTangents( mat, testF( 0.14 ), 3e-3, state, "cone, second increment" );
}

void testTangentAtTheApex()
{
  const auto props = properties( 5., 30., 20., 100., 0.05, 0.5 );
  Mat        mat( props.data(), props.size(), 1 );
  const auto state0 = freshState( mat );

  // a tension trial state beyond the apex, with a small deviatoric part
  const Tensor33d       F   = stretch( 1.01, 1.012, 1.0105 );
  const auto            res = evaluate( mat, F, 2e-3, state0 );
  const Eigen::Matrix3d tau = toEigen( res.tau );
  throwExceptionOnFailure( ( tau - tau.trace() / 3.0 * Eigen::Matrix3d::Identity() ).norm() <
                             1e-9 * std::abs( tau.trace() ),
                           "the tangent check must be at the apex" + where );
  checkTangents( mat, F, 2e-3, state0, "apex" );
}

void testNearTheVertex()
{
  // tension with shear, just below the apex pressure (from an MPM run, perfectly plastic cohesion): the solution
  // lies on the cone with a small deviator. The apex state here would violate its subdifferential condition, and a
  // return that falls back to it makes the response discontinuous in F.
  const std::array< double, 11 > props = { K, G, 5.0, 30.0, 10.0, 0.0, 0.005, 0.99, 2.0, 1.5, 1.0 };
  Mat                            mat( props.data(), props.size(), 1 );
  Tensor33d                      F = stretch( 1.0056742731167041, 0.99934189820798636, 1.0 );
  F( 0, 1 )                        = -0.0045365520453929751;
  F( 1, 0 )                        = -0.0055641206783383192;

  const auto            res = evaluate( mat, F, 0.0, freshState( mat ) );
  const Eigen::Matrix3d tau = toEigen( res.tau );
  const double          dev = ( tau - tau.trace() / 3.0 * Eigen::Matrix3d::Identity() ).norm();
  throwExceptionOnFailure( dev > 1e-2,
                           MakeString() << "the state must stay on the cone, |dev tau| = " << dev << where );
  throwExceptionOnFailure( std::abs( yieldFunction( res.tau, stateValue( mat, res.state, "alphaP" ), 5., 30., 0.0 ) ) <
                             1e-8,
                           "on the cone" + where );
  checkTangents( mat, F, 0.0, freshState( mat ), "near the vertex" );
}

void testLargeIncrement()
{
  // a far-off iterate of a global Newton scheme (from an MPM run): a plastic increment of order one, beyond the
  // range of the plain series of the tensor exponential. The return must still end on the cone.
  const auto props = properties( 5., 30., 10., 0.0 );
  Mat        mat( props.data(), props.size(), 1 );
  Tensor33d  F                                              = stretch( 1.6491543352545375, 0.67203941841161841, 1.0 );
  F( 0, 1 )                                                 = 0.0029042753684767816;
  F( 1, 0 )                                                 = 0.17121287455412429;
  std::vector< double > state                               = freshState( mat );
  double*               Fp                                  = mat.getStateView( "Fp", state.data() ).stateLocation;
  Fp[0]                                                     = 0.99901841034324146;
  Fp[4]                                                     = 1.001084449260059;
  Fp[8]                                                     = 1.0003438151045874;
  *mat.getStateView( "alphaP", state.data() ).stateLocation = 0.0025121126871932276;

  const auto   res    = evaluate( mat, F, 0.0, state );
  const double alphaP = stateValue( mat, res.state, "alphaP" );
  // reference: the previous, independent implementation (return map in the principal elastic log strains)
  throwExceptionOnFailure( std::abs( res.tau( 0, 0 ) + 66.70672101670 ) < 1e-8 &&
                             std::abs( res.tau( 1, 1 ) + 400.5796457121 ) < 1e-8 &&
                             std::abs( alphaP - 0.9634808934576 ) < 1e-11,
                           MakeString() << "large increment: tau_xx " << res.tau( 0, 0 ) << ", tau_yy "
                                        << res.tau( 1, 1 ) << ", alphaP " << alphaP << where );
  throwExceptionOnFailure( std::abs( yieldFunction( res.tau, alphaP, 5., 30., 0.0 ) ) < 1e-8,
                           "a large increment must still return to the cone" + where );
}

void testGradientEnhancedDamage()
{
  const double c0 = 5., phi = 30., psi = 20., H = 100., epsF = 0.05;
  const auto   props = properties( c0, phi, psi, H, epsF, 0.5 );
  Mat          mat( props.data(), props.size(), 1 );

  // the local variable grows with the dilatant plastic flow and is exported as the source L
  const auto res = evaluate( mat, testF( 0.1 ), 1e-2, freshState( mat ) );
  throwExceptionOnFailure( res.L > 0.0 && res.L == stateValue( mat, res.state, "alphaD" ),
                           "L must be the (positive) local damage variable" + where );

  const double omega = stateValue( mat, res.state, "omega" );
  const double aw    = 0.5 * 1e-2 + 0.5 * res.L;
  throwExceptionOnFailure( std::abs( omega - ( 1. - std::exp( -aw / epsF ) ) ) < 1e-14, "omega(alpha_w)" + where );

  // damage scales the effective stress
  const auto props0 = properties( c0, phi, psi, H, 1e10, 0.5 );
  Mat        undamaged( props0.data(), props0.size(), 1 );
  const auto res0 = evaluate( undamaged, testF( 0.1 ), 1e-2, freshState( undamaged ) );
  throwExceptionOnFailure( checkIfEqual( res.tau, Tensor33d( ( 1. - omega ) * res0.tau ), 1e-9 * norm( res0.tau ) ),
                           "tau != (1 - omega) tau_eff" + where );

  // the nonlocal field drives the damage ...
  throwExceptionOnFailure( norm( res.t.dTau_dN ) > 0.0, "tau must depend on the nonlocal field" + where );

  // ... irreversibly: a smaller nonlocal field in the next increment does not heal the material
  const auto later = evaluate( mat, testF( 0.1 ), 0.0, res.state );
  throwExceptionOnFailure( stateValue( mat, later.state, "omega" ) >= omega, "damage must not decrease" + where );
}

void testElasticUnloading()
{
  const auto props = properties( 5., 30., 10. );
  Mat        mat( props.data(), props.size(), 1 );

  const auto   loaded   = evaluate( mat, testF( 0.1 ), 0.0, freshState( mat ) );
  const double alphaP   = stateValue( mat, loaded.state, "alphaP" );
  const auto   unloaded = evaluate( mat, testF( 0.08 ), 0.0, loaded.state );

  throwExceptionOnFailure( stateValue( mat, unloaded.state, "alphaP" ) == alphaP,
                           "a small reverse increment must unload elastically" + where );
}

void testCumulativeDissipation()
{
  const auto props = properties( 5., 30., 10., 100., 0.05, 0.5 );
  Mat        mat( props.data(), props.size(), 1 );

  // the incoming dissipation is carried over and incremented by a plastic, damaging step ...
  const double carried = 3.0;
  const auto   loaded  = evaluate( mat, testF( 0.1 ), 1e-2, freshState( mat ), carried );
  throwExceptionOnFailure( stateValue( mat, loaded.state, "alphaP" ) > 0.0, "the step must yield" + where );
  throwExceptionOnFailure( loaded.dissipation > carried, "plastic flow and damage must dissipate" + where );

  // ... and left unchanged by an elastic one
  const auto elastic = evaluate( mat, testF( 0.08 ), 1e-2, loaded.state, loaded.dissipation );
  throwExceptionOnFailure( std::abs( elastic.dissipation - loaded.dissipation ) < 1e-14,
                           "an elastic step must not dissipate" + where );
}

void testFactoryAndValidation()
{
  const auto                                                    props = properties( 5., 30., 10. );
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > mat(
    MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::
      createMaterial( "GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER", props.data(), props.size(), 1 ) );
  throwExceptionOnFailure( mat != nullptr, "factory" + where );
  throwExceptionOnFailure( mat->getNumberOfRequiredStateVars() == 13,
                           "state layout: Fp, alphaP, alphaD, kappa, omega" + where );
  throwExceptionOnFailure( mat->getDensity( nullptr ) == props[10], "density" + where );

  bool threw = false;
  try {
    const auto bad = properties( 5., 20., 30. ); // dilatancy above friction
    Mat        m( bad.data(), bad.size(), 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "psi > phi must be rejected" + where );

  auto rejects = [&]( std::vector< double > p, const std::string& what ) {
    bool thrown = false;
    try {
      Mat m( p.data(), p.size(), 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, what + " must be rejected" + where );
  };
  const auto valid = properties( 5., 30., 10. );
  rejects( std::vector< double >( valid.begin(), valid.begin() + 9 ), "a too short property array" );
  for ( auto [idx, value, what] : std::vector< std::tuple< int, double, std::string > >{ { 0, 0.0, "K = 0" },
                                                                                         { 1, -1.0, "G < 0" },
                                                                                         { 7, 1.0, "maxDamage = 1" },
                                                                                         { 7, -0.1, "maxDamage < 0" },
                                                                                         { 2, 0.0, "c0 = 0" },
                                                                                         { 6, 0.0, "epsF = 0" } } ) {
    std::vector< double > p( valid.begin(), valid.end() );
    p[idx] = value;
    rejects( p, what );
  }
}

void testFailurePaths()
{
  // a degenerate deformation gradient cannot be decomposed
  {
    const auto props = properties( 5., 30., 10. );
    Mat        mat( props.data(), props.size(), 1 );
    Tensor33d  F = stretch( 1.0, 1.0, 1.0 );
    F( 2, 0 ) = F( 2, 1 ) = F( 2, 2 ) = 0.0;
    bool threw                        = false;
    try {
      evaluate( mat, F, 0.0, freshState( mat ) );
    }
    catch ( const Marmot::StressUpdateFailed& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "a singular F must raise StressUpdateFailed" + where );
  }

  // without dilatancy there is no apex to return to: hydrostatic tension beyond the cone has no admissible state
  {
    const auto props = properties( 5., 30., 0.0 );
    Mat        mat( props.data(), props.size(), 1 );
    bool       threw = false;
    try {
      evaluate( mat, stretch( 1.01, 1.01, 1.01 ), 0.0, freshState( mat ) );
    }
    catch ( const Marmot::StressUpdateFailed& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "a return map without solution must raise StressUpdateFailed" + where );
  }

  // near-apex tension with a small deviatoric part: the cone return would reverse the deviator, so the state must
  // end on the apex
  {
    const auto props = properties( 5., 30., 20., 100. );
    Mat        mat( props.data(), props.size(), 1 );
    const auto res = evaluate( mat, stretch( 1.01, 1.01, 1.0101 ), 0.0, freshState( mat ) );

    const Eigen::Matrix3d tau = toEigen( res.tau );
    const double          p   = tau.trace() / 3.0;
    throwExceptionOnFailure( ( tau - p * Eigen::Matrix3d::Identity() ).norm() < 1e-9 * std::abs( p ),
                             "a near-apex trial state must return to the apex" + where );
  }

  // the density is optional in the card, but asking for it without one is an error
  {
    const auto props = properties( 5., 30., 10. );
    Mat        mat( props.data(), 10, 1 );
    bool       threw = false;
    try {
      mat.getDensity( nullptr );
    }
    catch ( const std::runtime_error& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "a missing density must be reported" + where );
  }
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testElasticRangeIsCompressibleNeoHooke,
                                                               testConsistencyAndFlowRule,
                                                               testApexReturn,
                                                               testObjectivity,
                                                               testTangentInThePlasticBranch,
                                                               testTangentAtTheApex,
                                                               testLargeIncrement,
                                                               testNearTheVertex,
                                                               testGradientEnhancedDamage,
                                                               testElasticUnloading,
                                                               testFactoryAndValidation,
                                                               testFailurePaths };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
