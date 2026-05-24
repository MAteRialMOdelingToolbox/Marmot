#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot::Testing;
using namespace Marmot::Solvers;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// -----------------------------------------------------------------------
// Material property helpers
// -----------------------------------------------------------------------

// Layout: [hyperelasticBase, onlyShearCreep, K, G, nMaxwell, (gamma1, tau1, ...)]
// hyperelasticBase: 0=NeoHooke, 1=Yeoh, 2=MooneyRivlin, 3=PenceGouNeoHooke

static std::vector< double > getElasticNeoHookeProps()
{
  // NeoHooke, full creep flag off, K=3500, G=1500, nMaxwell=0
  return { 0.0, 0.0, 3500.0, 1500.0, 0.0 };
}

static std::vector< double > getViscoelasticNeoHookeProps()
{
  // NeoHooke, full creep, K=3500, G=1500, nMaxwell=1, gamma=0.3, tau=10
  return { 0.0, 0.0, 3500.0, 1500.0, 1.0, 0.3, 10.0 };
}

// Helper to create a single-step deformation solver
static MarmotMaterialPointSolverFiniteStrain makeSolver( const std::string& matName, std::vector< double >& matProps )
{
  auto solveropts = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  return MarmotMaterialPointSolverFiniteStrain( matName,
                                                matProps.data(),
                                                static_cast< int >( matProps.size() ),
                                                solveropts );
}

// Apply a single displacement-controlled step
static MarmotMaterialPointSolverFiniteStrain::Step makeStep( const Tensor33d& gradUIncrement,
                                                             double           timeStart,
                                                             double           timeEnd,
                                                             double           dT )
{
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = gradUIncrement;
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.timeStart                   = timeStart;
  step.timeEnd                     = timeEnd;
  step.dTStart                     = dT;
  step.dTMax                       = dT;
  return step;
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

// Test I-1: F=I gives zero Kirchhoff stress for elastic NeoHooke
void testUndeformedResponse()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  solver.addStep( makeStep( Tensor33d( 0.0 ), 0.0, 1.0, 1.0 ) );
  solver.solve();

  Tensor33d stressTarget( 0.0 );
  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, stressTarget, 1e-10 ),
                           "I-1: Undeformed configuration - stress should be zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Uniaxial stretch F_11=1.1 matches NeoHooke reference (elastic, nMaxwell=0)
void testUniaxialElasticResponse()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto finalStress = solver.getHistory().back().stress;

  // Reference: NeoHooke K=3500, G=1500, lambda=2500, F=diag(1.1,1,1)
  // PK2_11 = G*(1 - 1/C11) + lambda/2*(C11-1)/C11, C11=1.21
  // tau_11 = F11^2 * PK2_11 = 577.5, tau_22 = tau_33 = 262.5
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 577.5;
  stressTarget( 1, 1 ) = 262.5;
  stressTarget( 2, 2 ) = 262.5;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-4 ),
                           "I-2: Uniaxial elastic response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-3: Kirchhoff stress tensor is symmetric for arbitrary deformation
void testStressTensorSymmetry()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.01;
  gradU( 0, 1 ) = 0.06;
  gradU( 0, 2 ) = -0.03;
  gradU( 1, 0 ) = 0.06;
  gradU( 1, 1 ) = 0.02;
  gradU( 1, 2 ) = 0.04;
  gradU( 2, 0 ) = -0.03;
  gradU( 2, 1 ) = 0.04;
  gradU( 2, 2 ) = -0.05;
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto tau = solver.getHistory().back().stress;
  throwExceptionOnFailure( checkIfEqual( tau, Fastor::transpose( tau ), 1e-10 ),
                           "I-3: Stress symmetry failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-4: Pure rotation gives zero Kirchhoff stress
void testPureRotationZeroStress()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps();

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d F( 0.0 );
    F( 0, 0 ) = cos( phi );
    F( 0, 1 ) = -sin( phi );
    F( 1, 0 ) = sin( phi );
    F( 1, 1 ) = cos( phi );
    F( 2, 2 ) = 1.0;

    auto solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();

    throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, Tensor33d( 0.0 ), 1e-10 ),
                             "I-4: Pure rotation (phi=" + std::to_string( phi_deg ) + ") should give zero stress in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-5: Objectivity: tau(Q*F) = Q * tau(F) * Q^T
void testObjectivity()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps();

  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.01;
  F( 0, 1 ) = 0.06;
  F( 0, 2 ) = -0.03;
  F( 1, 0 ) = 0.06;
  F( 1, 1 ) = 1.02;
  F( 1, 2 ) = 0.04;
  F( 2, 0 ) = -0.03;
  F( 2, 1 ) = 0.04;
  F( 2, 2 ) = 0.95;

  auto solverRef = makeSolver( matName, matProps );
  solverRef.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
  solverRef.solve();
  Tensor33d tauRef = solverRef.getHistory().back().stress;

  for ( int phi_deg = 30; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.0;

    Tensor33d F_rot  = einsum< ik, kj, to_ij >( Q, F );
    auto      solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F_rot - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();
    Tensor33d tauRot = solver.getHistory().back().stress;

    // Expected: tau_rot = Q * tau_ref * Q^T
    Tensor33d tauExpected = einsum< iI, IJ, jJ, to_ij >( Q, tauRef, Q );

    throwExceptionOnFailure( checkIfEqual( tauRot, tauExpected, 1e-8 ),
                             "I-5: Objectivity failed for phi=" + std::to_string( phi_deg ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-6: Viscoelastic relaxation — long-term stress approaches (1-gamma)*sigma_elastic
void testViscoelasticRelaxation()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getViscoelasticNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  // Step 1: Apply instantaneous deformation (tiny time so almost no relaxation)
  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1e-4, 1e-4 ) );

  // Step 2: Hold deformation for 1000 >> tau=10 (fully relaxed)
  solver.addStep( makeStep( Tensor33d( 0.0 ), 1e-4, 1000.0, 100.0 ) );

  solver.solve();

  auto stressRelaxed = solver.getHistory().back().stress;

  // After full relaxation, tau = (1-gamma)*tau_elastic
  // NeoHooke K=3500, G=1500, F=diag(1.1,1,1): tau_11=577.5, tau_22=262.5
  const double gamma       = 0.3;
  const double longTermFac = 1.0 - gamma;

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = longTermFac * 577.5;
  stressTarget( 1, 1 ) = longTermFac * 262.5;
  stressTarget( 2, 2 ) = longTermFac * 262.5;

  throwExceptionOnFailure( checkIfEqual( stressRelaxed, stressTarget, 1.0 ),
                           "I-6: Viscoelastic relaxation - long-term stress wrong in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-7: Substepped variant gives the same result as the regular variant
void testSubsteppedConsistency()
{
  auto matPropsElas = getElasticNeoHookeProps();

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.05;
  gradU( 1, 1 ) = 0.02;
  gradU( 2, 2 ) = -0.01;

  auto step = makeStep( gradU, 0.0, 1.0, 1.0 );

  // Regular variant
  auto solverRegular = makeSolver( "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY", matPropsElas );
  solverRegular.addStep( step );
  solverRegular.solve();
  Tensor33d stressRegular = solverRegular.getHistory().back().stress;

  // Substepped variant
  std::vector< double > matPropsSub;
  matPropsSub.push_back( 1.0 ); // one substep
  for ( auto& prop : matPropsElas ) {
    matPropsSub.push_back( prop );
  }
  auto solverSub = makeSolver( "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED", matPropsSub );
  solverSub.addStep( step );
  solverSub.solve();
  Tensor33d stressSub = solverSub.getHistory().back().stress;

  throwExceptionOnFailure( checkIfEqual( stressRegular, stressSub, 1e-8 ),
                           "I-7: Substepped variant gives different result than regular variant in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,      // I-1: F=I gives zero stress
    testUniaxialElasticResponse, // I-2: Uniaxial elastic stretch reference values
    testStressTensorSymmetry,    // I-3: Kirchhoff stress symmetry
    testPureRotationZeroStress,  // I-4: Pure rotation gives zero stress
    testObjectivity,             // I-5: Objectivity tau(Q*F) = Q*tau(F)*Q^T
    testViscoelasticRelaxation,  // I-6: Viscoelastic relaxation
    testSubsteppedConsistency,   // I-7: Substepped == regular
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
