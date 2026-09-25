#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

using Mat = GradientEnhancedFiniteStrainDruckerPrager;

namespace {

  constexpr double K = 3500., G = 1500.;

  // K, G, c0, phi, psi, H, As, epsF, omegaMax, l, m, rho
  std::array< double, 12 > properties( double c0,
                                       double phi,
                                       double psi,
                                       double H    = 100.,
                                       double epsF = 1e10,
                                       double m    = 0.0 )
  {
    return { K, G, c0, phi, psi, H, 0.0, epsF, 0.99, 2.0, m, 2.4e-9 };
  }

  struct Result {
    Tensor33d                   tau;
    double                      L;
    std::vector< double >       state;
    Mat::AlgorithmicModuli< 3 > t;
  };

  // one call on a COPY of the given state (the history is not committed)
  Result evaluate( const Mat& mat, const Tensor33d& F, double N, const std::vector< double >& stateIn )
  {
    std::vector< double >          state = stateIn;
    Mat::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0., 0., 0., 0., state.data() );
    Mat::AlgorithmicModuli< 3 >    t;
    mat.computeStress( response, t, { F, N }, { 0.0, 1.0 } );
    return { response.tau, response.L, state, t };
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

void testTangentInThePlasticBranch()
{
  const auto props = properties( 5., 30., 10., 100., 0.05, 0.5 ); // plasticity AND damage active
  Mat        mat( props.data(), props.size(), 1 );

  const std::vector< double > state0 = freshState( mat );
  const Tensor33d             F      = testF( 0.1 );
  const double                N      = 2e-3;
  const auto                  res    = evaluate( mat, F, N, state0 );

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
          throwExceptionOnFailure( std::abs( res.t.dTau_dF( i, j, k, l ) - num ) < 1e-3 * ( 1. + std::abs( num ) ),
                                   "dTau_dF inconsistent" + where );
        }
      const double numL = ( rp.L - rm.L ) / ( 2 * h );
      throwExceptionOnFailure( std::abs( res.t.dL_dF( k, l ) - numL ) < 1e-4 * ( 1. + std::abs( numL ) ),
                               "dL_dF inconsistent" + where );
    }

  const double hN = 1e-7;
  const auto   rp = evaluate( mat, F, N + hN, state0 );
  const auto   rm = evaluate( mat, F, N - hN, state0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      const double num = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2 * hN );
      throwExceptionOnFailure( std::abs( res.t.dTau_dN( i, j ) - num ) < 1e-3 * ( 1. + std::abs( num ) ),
                               "dTau_dN inconsistent" + where );
    }
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

void testFactoryAndValidation()
{
  const auto                                                    props = properties( 5., 30., 10. );
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > mat(
    MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::
      createMaterial( "GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER", props.data(), props.size(), 1 ) );
  throwExceptionOnFailure( mat != nullptr, "factory" + where );
  throwExceptionOnFailure( mat->getNumberOfRequiredStateVars() == 13,
                           "state layout: Fp, alphaP, alphaD, kappa, omega" + where );
  throwExceptionOnFailure( mat->getDensity( nullptr ) == props[11], "density" + where );

  bool threw = false;
  try {
    const auto bad = properties( 5., 20., 30. ); // dilatancy above friction
    Mat        m( bad.data(), bad.size(), 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "psi > phi must be rejected" + where );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testElasticRangeIsCompressibleNeoHooke,
                                                               testConsistencyAndFlowRule,
                                                               testApexReturn,
                                                               testObjectivity,
                                                               testTangentInThePlasticBranch,
                                                               testGradientEnhancedDamage,
                                                               testElasticUnloading,
                                                               testFactoryAndValidation };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
