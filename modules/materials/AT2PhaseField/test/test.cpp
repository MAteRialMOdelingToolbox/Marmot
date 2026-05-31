#include "Marmot/AT2PhaseField.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Dense>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot;

// ─────────────────────────────────────────────────────────────────────────────
// Helper: create a fresh material + zero-initialised state variable vector.
// Material properties: E, nu, Gc, l
// ─────────────────────────────────────────────────────────────────────────────
static std::pair< AT2PhaseField, std::vector< double > > makeMaterial( const std::vector< double >& props )
{
  AT2PhaseField mat( props.data(), static_cast< int >( props.size() ), 1 );

  int                   nStateVars = mat.getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVars, 0.0 );
  return { std::move( mat ), std::move( stateVars ) };
}

// Convenience typedef
using Mat1 = MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >;
using Res1 = Mat1::response;
using Tan1 = Mat1::tangents;
using Inc1 = Mat1::increment;

// ─────────────────────────────────────────────────────────────────────────────
// Test 1 – Undamaged elastic response (phi = 0)
// Expected: stress = C:eps, KLocal = (2l/Gc)*H, c = l^2
// ─────────────────────────────────────────────────────────────────────────────
void testUndamagedElasticResponse()
{
  const double                E = 20000., nu = 0.25, Gc = 2.7, l = 0.01;
  const std::vector< double > properties = std::vector< double >{ E, nu, Gc, l };
  auto [mat, stateVars]                  = makeMaterial( properties );

  const Matrix6d C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

  // Apply a uniaxial strain increment, no phase-field
  Vector6d dEps      = Vector6d::Zero();
  dEps( 0 )          = 1e-3;
  const Vector6d eps = dEps; // starting from zero strain

  Res1 res;
  Tan1 tan;
  Inc1 inc;
  inc.dStrain = dEps;
  inc.K( 0 )  = 0.0; // phi = 0
  inc.dK( 0 ) = 0.0;
  inc.time    = 0.0;
  inc.dT      = 1.0;

  res.stress    = Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  mat.computeStress( res, tan, inc );

  // ── stress should equal C:eps ──
  const Vector6d stressExpected = C * eps;
  throwExceptionOnFailure( checkIfEqual< double >( res.stress, stressExpected, 1e-10 ),
                           "undamaged stress != C:eps in " + std::string( __PRETTY_FUNCTION__ ) );

  // ── gradient coefficient c = l^2 ──
  throwExceptionOnFailure( checkIfEqual( res.c( 0 ), l * l, 1e-15 ),
                           "c != l^2 in " + std::string( __PRETTY_FUNCTION__ ) );

  // ── KLocal = (2l/Gc)*(1-phi)*H  with phi=0 and H = psiPlus ──
  const double psiPlus      = 0.5 * eps.dot( C * eps );
  const double KLocalExpect = ( 2. * l / Gc ) * 1.0 * psiPlus;
  throwExceptionOnFailure( checkIfEqual( res.KLocal( 0 ), KLocalExpect, 1e-10 ),
                           "KLocal mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2 – Degraded stress with phi = 0.5
// Expected: stress = (1-phi)^2 * C:eps = 0.25 * C:eps
// ─────────────────────────────────────────────────────────────────────────────
void testDegradedStress()
{

  const double                E = 20000., nu = 0.25, Gc = 2.7, l = 0.01;
  const std::vector< double > properties = std::vector< double >{ E, nu, Gc, l };
  auto [mat, stateVars]                  = makeMaterial( properties );

  const Matrix6d C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

  Vector6d dEps = Vector6d::Zero();
  dEps( 0 )     = 1e-3;

  const double phi = 0.5;
  const double g   = ( 1. - phi ) * ( 1. - phi ); // = 0.25

  Res1 res;
  Tan1 tan;
  Inc1 inc;
  inc.dStrain = dEps;
  inc.K( 0 )  = phi;
  inc.dK( 0 ) = 0.0;
  inc.time    = 0.0;
  inc.dT      = 1.0;

  res.stress    = Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  mat.computeStress( res, tan, inc );

  const Vector6d stressExpected = g * ( C * dEps );
  throwExceptionOnFailure( checkIfEqual< double >( res.stress, stressExpected, 1e-10 ),
                           "degraded stress != g*C:eps in " + std::string( __PRETTY_FUNCTION__ ) );

  // KLocal should use (1-phi) factor
  const double psiPlus      = 0.5 * dEps.dot( C * dEps );
  const double KLocalExpect = ( 2. * l / Gc ) * ( 1. - phi ) * psiPlus;
  throwExceptionOnFailure( checkIfEqual( res.KLocal( 0 ), KLocalExpect, 1e-10 ),
                           "KLocal mismatch for phi=0.5 in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3 – Irreversibility of the crack-driving force history variable H
// Apply a large strain (step 1), then a smaller strain (step 2).
// H must not decrease.
// ─────────────────────────────────────────────────────────────────────────────
void testIrreversibility()
{
  const double                E = 20000., nu = 0.25, Gc = 2.7, l = 0.01;
  const std::vector< double > properties = std::vector< double >{ E, nu, Gc, l };
  auto [mat, stateVars]                  = makeMaterial( properties );

  // Step 1: large strain
  Vector6d dEps1 = Vector6d::Zero();
  dEps1( 0 )     = 2e-3;

  {
    Res1 res;
    Tan1 tan;
    Inc1 inc;
    inc.dStrain   = dEps1;
    inc.K( 0 )    = 0.0;
    inc.dK( 0 )   = 0.0;
    inc.time      = 0.0;
    inc.dT        = 1.0;
    res.stress    = Vector6d::Zero();
    res.KLocal    = Eigen::Vector< double, 1 >::Zero();
    res.c         = Eigen::Vector< double, 1 >::Zero();
    res.stateVars = stateVars.data();
    mat.computeStress( res, tan, inc );
  }

  // Capture H after step 1 from the state variable vector
  const double& H_after1 = stateVars[0]; // "maxCrackDrivingForce" is the first state var

  // Step 2: smaller strain increment that brings total strain back down
  Vector6d dEps2 = Vector6d::Zero();
  dEps2( 0 )     = -1e-3; // unloading: total strain becomes 1e-3

  {
    Res1 res;
    Tan1 tan;
    Inc1 inc;
    inc.dStrain   = dEps2;
    inc.K( 0 )    = 0.0;
    inc.dK( 0 )   = 0.0;
    inc.time      = 1.0;
    inc.dT        = 1.0;
    res.stress    = Vector6d::Zero();
    res.KLocal    = Eigen::Vector< double, 1 >::Zero();
    res.c         = Eigen::Vector< double, 1 >::Zero();
    res.stateVars = stateVars.data();
    mat.computeStress( res, tan, inc );
  }

  const double& H_after2 = stateVars[0];

  // H must not decrease
  throwExceptionOnFailure( H_after2 >= H_after1 - 1e-15,
                           "irreversibility violated: H decreased on unloading in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4 – Tangent consistency (numerical vs. analytical)
// Compare dStressddStrain, dStressddK, dKLocalddStrain, dKLocalddK
// using centred finite differences.
// ─────────────────────────────────────────────────────────────────────────────
void testTangentConsistency()
{
  const double                E = 20000., nu = 0.25, Gc = 2.7, l = 0.01;
  const std::vector< double > properties = std::vector< double >{ E, nu, Gc, l };
  auto [mat, stateVars]                  = makeMaterial( properties );

  const std::vector< double >
    stateVarsOld = stateVars; // capture original state variables to reset after each perturbation

  const double phi = 0.3;
  const double dK  = 0.0;

  // Base increment
  Vector6d dEps = Vector6d::Zero();
  dEps( 0 )     = 1e-3;
  dEps( 3 )     = 5e-4;

  // Compute base response and analytical tangent
  auto evalAt = [&]( const Vector6d& dE, double K_val, std::vector< double >& sv ) -> std::pair< Vector6d, double > {
    Res1 res;
    Tan1 tan;
    Inc1 inc;
    inc.dStrain   = dE;
    inc.K( 0 )    = K_val;
    inc.dK( 0 )   = dK;
    inc.time      = 0.0;
    inc.dT        = 1.0;
    res.stress    = Vector6d::Zero();
    res.KLocal    = Eigen::Vector< double, 1 >::Zero();
    res.c         = Eigen::Vector< double, 1 >::Zero();
    res.stateVars = sv.data();
    mat.computeStress( res, tan, inc );
    return { res.stress, res.KLocal( 0 ) };
  };

  // Get analytical tangent at base point
  Res1 resBase;
  Tan1 tanBase;
  Inc1 incBase;
  incBase.dStrain   = dEps;
  incBase.K( 0 )    = phi;
  incBase.dK( 0 )   = dK;
  incBase.time      = 0.0;
  incBase.dT        = 1.0;
  resBase.stress    = Vector6d::Zero();
  resBase.KLocal    = Eigen::Vector< double, 1 >::Zero();
  resBase.c         = Eigen::Vector< double, 1 >::Zero();
  resBase.stateVars = stateVars.data();
  mat.computeStress( resBase, tanBase, incBase );

  const double h = 1e-7; // finite-difference step

  // ── dStressddStrain ──
  for ( int j = 0; j < 6; j++ ) {
    Vector6d dEp = dEps, dEm = dEps;
    dEp( j ) += h;
    dEm( j ) -= h;
    stateVars     = stateVarsOld; // reset state variables before each perturbation to avoid accumulation of changes
    auto [sp, kp] = evalAt( dEp, phi, stateVars );
    stateVars     = stateVarsOld; // reset again before the second perturbation
    auto [sm, km] = evalAt( dEm, phi, stateVars );
    const Vector6d dSdE_num = ( sp - sm ) / ( 2. * h );
    const Vector6d dSdE_ana = tanBase.dStressddStrain.col( j );
    throwExceptionOnFailure( checkIfEqual< double >( dSdE_ana, dSdE_num, 1e-5 ),
                             "dStressddStrain col " + std::to_string( j ) + " mismatch in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // ── dStressddK (col 0) ──
  stateVars               = stateVarsOld; // reset state variables before perturbation
  auto [sp, kp]           = evalAt( dEps, phi + h, stateVars );
  stateVars               = stateVarsOld; // reset again before second perturbation
  auto [sm, km]           = evalAt( dEps, phi - h, stateVars );
  const Vector6d dSdK_num = ( sp - sm ) / ( 2. * h );
  const Vector6d dSdK_ana = tanBase.dStressddK.col( 0 );
  throwExceptionOnFailure( checkIfEqual< double >( dSdK_ana, dSdK_num, 1e-5 ),
                           "dStressddK mismatch in " + std::string( __PRETTY_FUNCTION__ ) );

  // ── dKLocalddStrain (row 0) ──
  for ( int j = 0; j < 6; j++ ) {
    Vector6d dEp = dEps, dEm = dEps;
    dEp( j ) += h;
    dEm( j ) -= h;
    stateVars             = stateVarsOld; // reset state variables before each perturbation
    auto [sp, kp]         = evalAt( dEp, phi, stateVars );
    stateVars             = stateVarsOld; // reset again before second perturbation
    auto [sm, km]         = evalAt( dEm, phi, stateVars );
    const double dKdE_num = ( kp - km ) / ( 2. * h );
    const double dKdE_ana = tanBase.dKLocalddStrain( 0, j );
    throwExceptionOnFailure( checkIfEqual( dKdE_ana, dKdE_num, 1e-5 ),
                             "dKLocalddStrain col " + std::to_string( j ) + " mismatch in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }

  // ── dKLocalddK (0,0) ──
  {
    stateVars             = stateVarsOld; // reset state variables before perturbation
    auto [sp, kp]         = evalAt( dEps, phi + h, stateVars );
    stateVars             = stateVarsOld; // reset again before second perturbation
    auto [sm, km]         = evalAt( dEps, phi - h, stateVars );
    const double dKdK_num = ( kp - km ) / ( 2. * h );
    const double dKdK_ana = tanBase.dKLocalddK( 0, 0 );
    throwExceptionOnFailure( checkIfEqual( dKdK_ana, dKdK_num, 1e-5 ),
                             "dKLocalddK mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}
// ─────────────────────────────────────────────────────────────────────────────
int main()
{
  const std::vector< std::function< void() > > tests = {
    testUndamagedElasticResponse,
    testDegradedStress,
    testIrreversibility,
    testTangentConsistency,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
