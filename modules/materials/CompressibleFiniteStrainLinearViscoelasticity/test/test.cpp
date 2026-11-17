#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::Math::Tensors::FastorTensors;
using namespace Marmot::Math::Tensors::FastorIndices;
using namespace Marmot::ContinuumMechanics;

// -----------------------------------------------------------------------
// Material property helpers
// -----------------------------------------------------------------------

// Layout: [hyperelasticBase, onlyShearCreep, <elastic params of the base potential...>, nMaxwell, (gamma1, tau1, ...)]
// hyperelasticBase: 0=NeoHooke (K, G), 1=Yeoh (C1, C2, C3, D1), 2=MooneyRivlin (C1, C2, D1),
// 3=PenceGouNeoHooke (K, G), 4=ArrudaBoyce (mu, lambdaL, D1)

// The elastic NeoHooke tests run for both onlyShearCreep=0 (full creep) and onlyShearCreep=1, since the flag
// selects different initial-stiffness paths in the constructor.
static std::vector< double > getElasticNeoHookeProps( double onlyShearCreep )
{
  // NeoHooke, K=3500, G=1500, nMaxwell=0
  return { 0.0, onlyShearCreep, 3500.0, 1500.0, 0.0 };
}

static std::vector< double > getViscoelasticNeoHookeFullCreepProps()
{
  // NeoHooke, full creep (onlyShearCreep off), K=3500, G=1500, nMaxwell=1, gamma=0.3, tau=10
  return { 0.0, 0.0, 3500.0, 1500.0, 1.0, 0.3, 10.0 };
}

static std::vector< double > getViscoelasticNeoHookeProps()
{
  // NeoHooke, onlyShearCreep on, K=3500, G=1500, nMaxwell=1, gamma=0.3, tau=10
  return { 0.0, 1.0, 3500.0, 1500.0, 1.0, 0.3, 10.0 };
}

// Layout: [hyperelasticBase=4 (ArrudaBoyce), onlyShearCreep, mu, lambdaL, D1, nMaxwell]
static std::vector< double > getElasticArrudaBoyceProps()
{
  // ArrudaBoyce, full creep flag off, mu=1500, lambdaL=3.0, D1=2/3500 (bulk modulus 3500), nMaxwell=0
  return { 4.0, 0.0, 1500.0, 3.0, 2.0 / 3500.0, 0.0 };
}

// Layout: [hyperelasticBase=4 (ArrudaBoyce), onlyShearCreep, mu, lambdaL, D1, nMaxwell, gamma1, tau1]
static std::vector< double > getViscoelasticArrudaBoyceShearCreepProps()
{
  // ArrudaBoyce, onlyShearCreep on, mu=1500, lambdaL=3.0, D1=2/3500, nMaxwell=1, gamma=0.3, tau=10
  return { 4.0, 1.0, 1500.0, 3.0, 2.0 / 3500.0, 1.0, 0.3, 10.0 };
}

static std::vector< double > getElasticPenceGouNeoHookeProps()
{
  // PenceGouNeoHooke, full creep flag off, K=3500, G=1500, nMaxwell=0
  return { 3.0, 0.0, 3500.0, 1500.0, 0.0 };
}

// Helper to create a single-step deformation solver
static Marmot::Materials::MaterialPointSolverFiniteStrain makeSolver( const std::string&     matName,
                                                                      std::vector< double >& matProps )
{
  auto solveropts = Marmot::Materials::MaterialPointSolverFiniteStrain::SolverOptions();
  return Marmot::Materials::MaterialPointSolverFiniteStrain( matName,
                                                             matProps.data(),
                                                             static_cast< int >( matProps.size() ),
                                                             solveropts );
}

// Apply a single displacement-controlled step
static Marmot::Materials::MaterialPointSolverFiniteStrain::Step makeStep( const Tensor33d& gradUIncrement,
                                                                          double           timeStart,
                                                                          double           timeEnd,
                                                                          double           dT )
{
  Marmot::Materials::MaterialPointSolverFiniteStrain::Step step;
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
void testUndeformedResponse( double onlyShearCreep )
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps( onlyShearCreep );
  auto              solver   = makeSolver( matName, matProps );

  solver.addStep( makeStep( Tensor33d( 0.0 ), 0.0, 1.0, 1.0 ) );
  solver.solve();

  Tensor33d stressTarget( 0.0 );
  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, stressTarget, 1e-10 ),
                           "I-1: Undeformed configuration - stress should be zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Uniaxial stretch F_11=1.1 matches NeoHooke reference (elastic, nMaxwell=0)
void testUniaxialElasticResponse( double onlyShearCreep )
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps( onlyShearCreep );
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto finalStress = solver.getHistory().back().stress;

  // Reference: NeoHooke K=3500, G=1500, F=diag(1.1,1,1): tau = G (b - I) + lambda/2 (J^2 - 1) I
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 577.5;
  stressTarget( 1, 1 ) = 262.5;
  stressTarget( 2, 2 ) = 262.5;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-4 ),
                           "I-2: Uniaxial elastic response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-3: Kirchhoff stress tensor is symmetric for arbitrary deformation
void testStressTensorSymmetry( double onlyShearCreep )
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps( onlyShearCreep );
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
void testPureRotationZeroStress( double onlyShearCreep )
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps( onlyShearCreep );

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
void testObjectivity( double onlyShearCreep )
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticNeoHookeProps( onlyShearCreep );

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

  // onlyShearCreep=1: only the deviatoric part of PK2 relaxes, PK2vol stays at
  // its instantaneous value -- tau_final = F (PK2vol + (1-gamma) PK2dev) F^T with gamma=0.3,
  // for the NeoHooke K=3500, G=1500, F=diag(1.1,1,1) response of testUniaxialElasticResponse.
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 525.525;
  stressTarget( 1, 1 ) = 283.977272727273;
  stressTarget( 2, 2 ) = 283.977272727273;

  throwExceptionOnFailure( checkIfEqual( stressRelaxed, stressTarget, 1.0 ),
                           "I-6: Viscoelastic relaxation - long-term stress wrong in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-7: Substepped variant gives the same result as the regular variant
void testSubsteppedConsistency( double onlyShearCreep )
{
  auto matPropsElas = getElasticNeoHookeProps( onlyShearCreep );

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

// Test I-8: F=I gives zero Kirchhoff stress for elastic ArrudaBoyce
void testArrudaBoyceUndeformedResponse()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticArrudaBoyceProps();
  auto              solver   = makeSolver( matName, matProps );

  solver.addStep( makeStep( Tensor33d( 0.0 ), 0.0, 1.0, 1.0 ) );
  solver.solve();

  Tensor33d stressTarget( 0.0 );
  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, stressTarget, 1e-10 ),
                           "I-8: ArrudaBoyce undeformed configuration - stress should be zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-9: ArrudaBoyce response matches the shared core-header potential directly,
// end-to-end through this material's autodiff::dual3rd-based
// computeEnergyDensityAndDerivatives wiring and its shared
// D1-based volumetric term (shared with Yeoh and MooneyRivlin) -- NOT a check of the potential
// formula itself (that is covered by
// TestMarmotEnergyDensityFunctions), but of the switch-case dispatch and
// derivative pipeline added for this base in this material.
void testArrudaBoyceMatchesSharedPotential()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticArrudaBoyceProps();
  const double      mu = matProps[2], lambdaL = matProps[3], D1 = matProps[4];
  auto              solver = makeSolver( matName, matProps );

  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.08;
  F( 1, 1 ) -= 0.02;
  F( 2, 2 ) -= 0.015;
  F( 0, 1 ) = 0.01;

  solver.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
  solver.solve();
  Tensor33d tau = solver.getHistory().back().stress;

  const Tensor33d C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );
  double          psiIso;
  Tensor33d       dPsiIso_dC;
  std::tie( psiIso, dPsiIso_dC ) = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric::
    FirstOrder::ArrudaBoyce8ChainPotential< double >( C, mu, lambdaL );

  const double    detC       = Fastor::determinant( C );
  const Tensor33d CInv       = Fastor::inverse( C );
  const Tensor33d dPsiVol_dC = ( ( detC - 1. ) / ( 2. * D1 ) ) * CInv;

  const Tensor33d PK2       = 2. * ( dPsiIso_dC + dPsiVol_dC );
  const Tensor33d tauTarget = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

  throwExceptionOnFailure( checkIfEqual( tau, tauTarget, 1e-6 ),
                           "I-9: ArrudaBoyce stress does not match the shared potential's directly-assembled "
                           "target in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-10: F=I gives zero Kirchhoff stress for elastic PenceGouNeoHooke (exercises the hyperelasticBase branches
// of computeEnergyDensityAndDerivatives other than the default NeoHooke one, e.g. the third-derivative computation
// via Marmot::Math::Differentiation::Automatic::ThirdOrder::d3f_dT3)
void testUndeformedResponsePenceGouNeoHooke()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticPenceGouNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  solver.addStep( makeStep( Tensor33d( 0.0 ), 0.0, 1.0, 1.0 ) );
  solver.solve();

  Tensor33d stressTarget( 0.0 );
  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, stressTarget, 1e-10 ),
                           "I-10: Undeformed configuration (PenceGouNeoHooke) - stress should be zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-11: Triaxial stretch for PenceGouNeoHooke matches an independently computed reference (PK2 = 2*dPsi/dC via
// Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::PenceGouPotentialB, converted to
// Kirchhoff stress)
void testTriaxialResponsePenceGouNeoHookeMatchesAnalytical()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getElasticPenceGouNeoHookeProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1;
  gradU( 1, 1 ) = -0.05;
  gradU( 2, 2 ) = 0.02;
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto finalStress = solver.getHistory().back().stress;

  const double K = matProps[2];
  const double G = matProps[3];
  Tensor33d    F = gradU + Spatial3D::I;
  Tensor33d    C = Marmot::ContinuumMechanics::Kinematics::DeformationMeasures::rightCauchyGreen( F );

  auto [psi,
        dPsi_dC] = Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::FirstOrder::PenceGouPotentialB( C,
                                                                                                                   K,
                                                                                                                   G );
  Tensor33d PK2  = 2. * dPsi_dC;
  Tensor33d stressTarget = Marmot::ContinuumMechanics::StressMeasures::KirchhoffStressFromPK2( PK2, F );

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-4 ),
                           "I-11: Triaxial PenceGouNeoHooke response doesn't match the analytical reference in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-12: Arruda-Boyce viscoelastic shear creep with onlyShearCreep=1 exercises the
// closed-form infinitesimal shear modulus branch G0 = mu/3 + 2*mu/(3*(1-1/lambdaL^2)) in
// initialCompliance = 0.5/G0 * I4.
// For a small simple-shear step (gamma12 << 1):
// 1) The instantaneous response matches the purely elastic Arruda-Boyce response and the
//    small-shear analytic shear stress G0 * gamma12, which directly validates that
//    initialCompliance : initialTangent = I4 (an erroneous G0=mu would produce a ~2.5% discrepancy).
// 2) The fully relaxed hold (t = 1000 >> tau) relaxes to (1 - gamma) * G0 * gamma12.
void testArrudaBoyceViscoelasticShearCreep()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getViscoelasticArrudaBoyceShearCreepProps();
  const double      mu       = matProps[2];
  const double      lambdaL  = matProps[3];
  const double      gamma    = matProps[6];
  const double      tau      = matProps[7];
  auto              solver   = makeSolver( matName, matProps );

  const double w0 = 1.0 - 1.0 / ( lambdaL * lambdaL );
  const double G0 = mu / 3.0 + 2.0 * mu / ( 3.0 * w0 );

  const double gamma12 = 1e-4;
  Tensor33d    gradU( 0.0 );
  gradU( 0, 1 ) = gamma12;

  // Step 1: Instantaneous loading step (dt << tau)
  solver.addStep( makeStep( gradU, 0.0, 1e-6, 1e-6 ) );

  // Step 2: Fully relaxed hold (t = 100 * tau >> tau)
  solver.addStep( makeStep( Tensor33d( 0.0 ), 1e-6, 100.0 * tau, 10.0 * tau ) );

  solver.solve();

  const auto& history = solver.getHistory();
  Tensor33d   tauInstant( 0.0 );
  for ( const auto& h : history ) {
    if ( std::abs( h.time - 1e-6 ) < 1e-12 ) {
      tauInstant = h.stress;
      break;
    }
  }
  const Tensor33d tauRelaxed = history.back().stress;

  const double targetInstant = G0 * gamma12;
  const double targetRelaxed = ( 1.0 - gamma ) * G0 * gamma12;

  // 1. Instantaneous shear stress must match analytic G0 * gamma12
  // (guards against initialCompliance using mu instead of G0, which would produce ~2.5% discrepancy)
  throwExceptionOnFailure( std::abs( tauInstant( 0, 1 ) - targetInstant ) < 1e-6,
                           "I-12: ArrudaBoyce shear-creep instantaneous shear stress does not match analytic G0 target "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // 2. Relaxed shear stress must match (1 - gamma) * G0 * gamma12
  throwExceptionOnFailure( std::abs( tauRelaxed( 0, 1 ) - targetRelaxed ) < 1e-6,
                           "I-12: ArrudaBoyce shear-creep relaxed shear stress does not match target in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-13: full creep (onlyShearCreep=0) with NeoHooke -- the whole PK2 stress relaxes, so the long-term
// stress is (1-gamma) times the instantaneous elastic response of testUniaxialElasticResponse. Exercises the
// constructor's inversion of the autodiff-based initial stiffness, which requires minor-symmetric tangents.
void testFullCreepViscoelasticRelaxation()
{
  const std::string matName  = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  auto              matProps = getViscoelasticNeoHookeFullCreepProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1e-4, 1e-4 ) );
  solver.addStep( makeStep( Tensor33d( 0.0 ), 1e-4, 1000.0, 100.0 ) );
  solver.solve();

  const double gamma = 0.3;
  Tensor33d    stressTarget( 0.0 );
  stressTarget( 0, 0 ) = ( 1.0 - gamma ) * 577.5;
  stressTarget( 1, 1 ) = ( 1.0 - gamma ) * 262.5;
  stressTarget( 2, 2 ) = ( 1.0 - gamma ) * 262.5;

  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress, stressTarget, 1.0 ),
                           "I-13: Full-creep viscoelastic relaxation - long-term stress wrong in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-14: for every hyperelastic base the analytical initial shear modulus G0 defines the deviatoric
// part of the initial stiffness dPK2/dE = 4 d2Psi/dCdC at C = I: Dev : dPK2/dE : Dev = 2 G0 Dev. This is
// what justifies initialCompliance = 0.5/G0 * I4 for onlyShearCreep = 1.
void testAnalyticalInitialShearModulusMatchesInitialStiffness()
{
  struct Exposed : CompressibleFiniteStrainLinearViscoelasticity {
    using CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity;
    using CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives;
    using CompressibleFiniteStrainLinearViscoelasticity::initialShearModulus;
  };

  // [hyperelasticBase, onlyShearCreep, <elastic params>, nMaxwell], expected G0
  const std::vector< std::pair< std::vector< double >, double > > cases = {
    { { 0.0, 1.0, 3500.0, 1500.0, 0.0 }, 1500.0 },                                   // NeoHooke: G
    { { 1.0, 1.0, 700.0, -50.0, 20.0, 2.0 / 3500.0, 0.0 }, 1400.0 },                 // Yeoh: 2 C1
    { { 2.0, 1.0, 500.0, 200.0, 2.0 / 3500.0, 0.0 }, 1400.0 },                       // MooneyRivlin: 2 (C1 + C2)
    { { 3.0, 1.0, 3500.0, 1500.0, 0.0 }, 1500.0 },                                   // PenceGouNeoHooke: G
    { { 4.0, 1.0, 1500.0, 3.0, 2.0 / 3500.0, 0.0 }, 1500.0 * 26.0 / ( 3.0 * 8.0 ) }, // ArrudaBoyce, lambdaL = 3
  };

  const Tensor3333d& Dev = Spatial3D::DeviatoricSymmetric;

  for ( const auto& [matProps, G0Expected] : cases ) {
    const Exposed material( matProps.data(), static_cast< int >( matProps.size() ), 0 );
    const double  G0 = material.initialShearModulus();

    throwExceptionOnFailure( checkIfEqual( G0, G0Expected, 1e-10 * G0Expected ),
                             "I-14: unexpected analytical initial shear modulus for base " +
                               std::to_string( static_cast< int >( matProps[0] ) ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );

    const Tensor3333d initialStiffness = 4. * std::get< 2 >( material.computeEnergyDensityAndDerivatives(
                                                Fastor::evaluate( Spatial3D::I ) ) );
    const Tensor3333d
      deviatoricStiffness = Fastor::einsum< ijmn, mnKL >( Dev, Fastor::einsum< ijmn, mnKL >( initialStiffness, Dev ) );
    const Tensor3333d difference = deviatoricStiffness - 2. * G0 * Dev;

    throwExceptionOnFailure( Fastor::max( Fastor::abs( difference ) ) < 1e-9 * G0,
                             "I-14: deviatoric initial stiffness is not 2 G0 Dev for base " +
                               std::to_string( static_cast< int >( matProps[0] ) ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// runs a NeoHooke test for full creep (onlyShearCreep=0) and shear creep only (onlyShearCreep=1)
static std::function< void() > forBothCreepModes( void ( *test )( double ) )
{
  return [test]() {
    for ( const double onlyShearCreep : { 0.0, 1.0 } )
      test( onlyShearCreep );
  };
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    forBothCreepModes( testUndeformedResponse ),              // I-1: F=I gives zero stress
    forBothCreepModes( testUniaxialElasticResponse ),         // I-2: Uniaxial elastic stretch reference values
    forBothCreepModes( testStressTensorSymmetry ),            // I-3: Kirchhoff stress symmetry
    forBothCreepModes( testPureRotationZeroStress ),          // I-4: Pure rotation gives zero stress
    forBothCreepModes( testObjectivity ),                     // I-5: Objectivity tau(Q*F) = Q*tau(F)*Q^T
    testViscoelasticRelaxation,                               // I-6: Viscoelastic relaxation
    forBothCreepModes( testSubsteppedConsistency ),           // I-7: Substepped == regular
    testArrudaBoyceUndeformedResponse,                        // I-8: ArrudaBoyce F=I gives zero stress
    testArrudaBoyceMatchesSharedPotential,                    // I-9: ArrudaBoyce matches shared potential
    testUndeformedResponsePenceGouNeoHooke,                   // I-10: F=I gives zero stress (PenceGouNeoHooke)
    testTriaxialResponsePenceGouNeoHookeMatchesAnalytical,    // I-11: PenceGouNeoHooke matches analytical reference
    testArrudaBoyceViscoelasticShearCreep,                    // I-12: ArrudaBoyce viscoelastic shear creep
    testFullCreepViscoelasticRelaxation,                      // I-13: full creep (onlyShearCreep=0) relaxation
    testAnalyticalInitialShearModulusMatchesInitialStiffness, // I-14: analytical G0 for every base
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
