#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

using Mat = GradientEnhancedCompressibleNeoHookeDamage;

namespace {

  // K, G, kappa0, kappaF, l, rho
  const std::array< double, 6 > props = { 3500., 1500., 1e-3, 1e-2, 2.5, 2.4e-9 };

  const Mat& material()
  {
    static const Mat mat( props.data(), props.size(), 1 );
    return mat;
  }

  struct Result {
    Tensor33d                   tau;
    double                      L;
    double                      kappa;
    Mat::AlgorithmicModuli< 3 > t;
  };

  // evaluate with a COPY of the history, so that repeated calls see the same state
  Result evaluate( const Tensor33d& F, double N, double kappaIn )
  {
    double                         kappa = kappaIn;
    Mat::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0., 0., 0., 0., &kappa );
    Mat::AlgorithmicModuli< 3 >    t;
    material().computeStress( response, t, { F, N }, { 0.0, 1.0 } );
    return { response.tau, response.L, kappa, t };
  }

  Tensor33d testF()
  {
    Tensor33d F = { { 1.012, 0.004, -0.002 }, { 0.003, 0.995, 0.006 }, { -0.001, 0.002, 1.004 } };
    return F;
  }

  const std::string where = " in " + std::string( __FILE__ );

} // namespace

void testNoDamageBelowThresholdIsNeoHooke()
{
  CompressibleNeoHooke                            reference( props.data(), 2, 1 );
  CompressibleNeoHooke::ConstitutiveResponse< 3 > r( Tensor33d( 0.0 ), 0., 0., nullptr );
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    t;
  reference.computeStress( r, t, { testF() }, { 0.0, 1.0 } );

  const auto res = evaluate( testF(), 0.5 * props[2], 0.0 );

  throwExceptionOnFailure( checkIfEqual( res.tau, r.tau, 1e-10 ),
                           "tau below threshold != CompressibleNeoHooke" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dTau_dF, t.dTau_dF, 1e-8 ),
                           "dTau_dF below threshold != CompressibleNeoHooke" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dTau_dN, Tensor33d( 0.0 ), 1e-14 ), "dTau_dN must vanish" + where );
}

void testUndeformedState()
{
  Tensor33d I;
  I.eye();
  const auto res = evaluate( I, 0.0, 0.0 );
  throwExceptionOnFailure( checkIfEqual( res.tau, Tensor33d( 0.0 ), 1e-12 ), "tau(I) must vanish" + where );
  throwExceptionOnFailure( checkIfEqual( res.L, 0.0, 1e-14 ), "L(I) must vanish" + where );
}

void testEquivalentStrainUniaxial()
{
  // small uniaxial stretch with lateral contraction nu: L -> axial strain
  const double eps = 1e-5;
  const double K = props[0], G = props[1];
  const double nu = ( 3. * K - 2. * G ) / ( 2. * ( 3. * K + G ) );
  Tensor33d    F( 0.0 );
  F( 0, 0 )      = 1. + eps;
  F( 1, 1 )      = 1. - nu * eps;
  F( 2, 2 )      = 1. - nu * eps;
  const auto res = evaluate( F, 0.0, 0.0 );
  throwExceptionOnFailure( std::abs( res.L - eps ) < 1e-3 * eps,
                           "L != uniaxial strain in the small-strain limit" + where );
}

void testConsistentTangentsLoading()
{
  const Tensor33d F     = testF();
  const double    N     = 3e-3; // > kappa0, > history -> loading branch
  const double    kappa = 2e-3;
  const auto      res   = evaluate( F, N, kappa );

  throwExceptionOnFailure( res.kappa == N, "history not updated on loading" + where );

  const double h = 1e-7;

  Tensor3333d dTau_dF_num( 0.0 );
  Tensor33d   dL_dF_num( 0.0 );
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Tensor33d Fp = F, Fm = F;
      Fp( k, l ) += h;
      Fm( k, l ) -= h;
      const auto rp = evaluate( Fp, N, kappa );
      const auto rm = evaluate( Fm, N, kappa );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          dTau_dF_num( i, j, k, l ) = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2 * h );
      dL_dF_num( k, l ) = ( rp.L - rm.L ) / ( 2 * h );
    }

  const double    hN          = 1e-9;
  const auto      rpN         = evaluate( F, N + hN, kappa );
  const auto      rmN         = evaluate( F, N - hN, kappa );
  const Tensor33d dTau_dN_num = Tensor33d( ( rpN.tau - rmN.tau ) / ( 2 * hN ) );
  const double    dL_dN_num   = ( rpN.L - rmN.L ) / ( 2 * hN );

  throwExceptionOnFailure( checkIfEqual( res.t.dTau_dF, dTau_dF_num, 1e-4 ), "dTau_dF inconsistent" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dTau_dN, dTau_dN_num, 1e-2 ), "dTau_dN inconsistent" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dL_dF, dL_dF_num, 1e-7 ), "dL_dF inconsistent" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dL_dN, dL_dN_num, 1e-10 ), "dL_dN inconsistent" + where );
}

void testUnloadingKeepsDamage()
{
  const double kappa = 5e-3;
  const auto   res   = evaluate( testF(), 2e-3, kappa ); // N below history -> elastic unloading
  const auto [D, _]  = material().damage( kappa );

  throwExceptionOnFailure( res.kappa == kappa, "history must not decrease" + where );
  throwExceptionOnFailure( checkIfEqual( res.t.dTau_dN, Tensor33d( 0.0 ), 1e-14 ),
                           "dTau_dN must vanish when unloading" + where );

  const auto undamaged = evaluate( testF(), 0.0, 0.0 );
  throwExceptionOnFailure( checkIfEqual( res.tau, Tensor33d( ( 1. - D ) * undamaged.tau ), 1e-10 ),
                           "tau != (1-D) tau0" + where );
}

void testObjectivity()
{
  const auto ref = evaluate( testF(), 4e-3, 0.0 );
  for ( int deg = 0; deg <= 180; deg += 45 ) {
    const double phi = Marmot::Math::degToRad( deg );
    Tensor33d    Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.;

    const Tensor33d QF  = einsum< ik, kj, to_ij >( Q, testF() );
    const auto      rot = evaluate( QF, 4e-3, 0.0 );

    const Tensor33d tauRotated = einsum< iI, IJ, jJ, to_ij >( Q, ref.tau, Q );
    throwExceptionOnFailure( checkIfEqual( rot.tau, tauRotated, 1e-9 ), "objectivity of tau failed" + where );
    throwExceptionOnFailure( checkIfEqual( rot.L, ref.L, 1e-12 ), "L must be invariant under rotations" + where );
  }
}

void testFactoryAndProperties()
{
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > mat(
    MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::
      createMaterial( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE", props.data(), props.size(), 7 ) );
  throwExceptionOnFailure( mat != nullptr, "factory did not create the material" + where );
  throwExceptionOnFailure( mat->getNumberOfRequiredStateVars() == 1, "expected one state variable (kappa)" + where );
  throwExceptionOnFailure( checkIfEqual( mat->getDensity( nullptr ), props[5] ), "density" + where );

  bool threw = false;
  try {
    std::array< double, 5 > bad = { 3500., 1500., 1e-2, 1e-3, 1. }; // kappaF < kappa0
    Mat                     m( bad.data(), bad.size(), 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "kappaF <= kappa0 must be rejected" + where );
}

void testCumulativeDissipation()
{
  // carried over, incremented by a damaging step, unchanged by an undamaged one
  const double                   carried = 2.0;
  double                         kappa   = 2e-3;
  Mat::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0., 0., 0., carried, &kappa );
  Mat::AlgorithmicModuli< 3 >    t;
  material().computeStress( response, t, { testF(), 3e-3 }, { 0.0, 1.0 } );
  throwExceptionOnFailure( response.dissipation > carried, "damage growth must dissipate" + where );

  const double afterDamage = response.dissipation;
  material().computeStress( response, t, { testF(), 1e-3 }, { 0.0, 1.0 } ); // below the history: no new damage
  throwExceptionOnFailure( std::abs( response.dissipation - afterDamage ) < 1e-14,
                           "a step without damage growth must not dissipate" + where );
}

void testPropertyValidation()
{
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
  rejects( { 3500., 1500., 1e-3, 1e-2 }, "a too short property array" );
  rejects( { 0.0, 1500., 1e-3, 1e-2, 1. }, "K = 0" );
  rejects( { 3500., -1., 1e-3, 1e-2, 1. }, "G < 0" );
  rejects( { 3500., 1500., 0.0, 1e-2, 1. }, "kappa0 = 0" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testNoDamageBelowThresholdIsNeoHooke,
                                                               testUndeformedState,
                                                               testEquivalentStrainUniaxial,
                                                               testConsistentTangentsLoading,
                                                               testUnloadingKeepsDamage,
                                                               testObjectivity,
                                                               testFactoryAndProperties,
                                                               testCumulativeDissipation,
                                                               testPropertyValidation };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
