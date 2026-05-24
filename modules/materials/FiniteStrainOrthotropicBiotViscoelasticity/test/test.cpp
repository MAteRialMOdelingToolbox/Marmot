#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot::Testing;
using namespace Marmot::Solvers;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// Material properties for an isotropic material (all directions equal)
// Layout: [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, nMaxwell]
// nMaxwell=0 means purely elastic, no Maxwell elements
static std::vector< double > getIsotropicProperties()
{
  const double E   = 20000.0;
  const double nu  = 0.25;
  const double G   = E / ( 2.0 * ( 1.0 + nu ) );
  const double nMW = 0;
  return { E, E, E, nu, nu, nu, G, G, G, nMW };
}

// Material properties for an orthotropic material
// Layout: [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, nMaxwell]
// nMaxwell=0 means purely elastic, no Maxwell elements
static std::vector< double > getOrthotropicProperties()
{
  return { 20000.0, 10000.0, 15000.0, 0.25, 0.35, 0.30, 4000.0, 6000.0, 5000.0, 0 };
}

// Material properties for an isotropic material with one Maxwell element for viscoelastic tests
// Layout: [E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, nMaxwell=1, gamma1, tau1]
static std::vector< double > getIsotropicViscoelasticProperties()
{
  const double E   = 20000.0;
  const double nu  = 0.25;
  const double G   = E / ( 2.0 * ( 1.0 + nu ) );
  const double nMW = 1;
  // Maxwell element: gamma=0.3 (relaxation amplitude), tau=10 (relaxation time)
  return { E, E, E, nu, nu, nu, G, G, G, nMW, 0.3, 10.0 };
}

// Helper to create a solver with the given material name and properties
static MarmotMaterialPointSolverFiniteStrain makeSolver( const std::string&     matName,
                                                         std::vector< double >& matProps,
                                                         const Tensor33d&       initialStress    = Tensor33d( 0.0 ),
                                                         const Eigen::VectorXd& initialStateVars = Eigen::VectorXd() )
{
  auto solveropts = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto solver     = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       matProps.data(),
                                                       static_cast< int >( matProps.size() ),
                                                       solveropts );
  return solver;
}

// Test I-1: Undeformed configuration F=I should yield zero Kirchhoff stress
void testUndeformedResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  auto solver = makeSolver( matName, matProps );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = Tensor33d( 0.0 );
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.timeStart                   = 0.0;
  step.timeEnd                     = 1.0;
  step.dTStart                     = 1.0;
  step.dTMax                       = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  Tensor33d stressTarget( 0.0 );

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-10 ),
                           "I-1: Undeformed configuration - Kirchhoff stress should be zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Small uniaxial stretch in direction 1 (isotropic, purely elastic)
// Reference values computed analytically from the Biot hyperelastic formulation:
// For F = diag(1+eps, 1, 1) with isotropic E=20000, nu=0.25:
//   tau_11 = (1+eps) * C_11 * eps, tau_22 = tau_33 = C_12 * eps
//   where C_11 = E*(1-nu)/((1+nu)*(1-2*nu)) = 24000
//         C_12 = E*nu/((1+nu)*(1-2*nu))     = 8000
void testSmallStrainUniaxialResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  auto solver = makeSolver( matName, matProps );

  const double eps = 0.001;

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = eps;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 1.0;
  step.dTMax                        = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  // tau_11 = (1+eps) * 24000 * eps = 1.001 * 24 = 24.024
  // tau_22 = tau_33 = 8000 * eps = 8.0
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = ( 1.0 + eps ) * 24000.0 * eps;
  stressTarget( 1, 1 ) = 8000.0 * eps;
  stressTarget( 2, 2 ) = 8000.0 * eps;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-8 ),
                           "I-2: Small strain uniaxial response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-3: Finite strain uniaxial stretch (isotropic, purely elastic)
// For F = diag(1.1, 1, 1) with isotropic E=20000, nu=0.25:
//   tau_11 = 2640, tau_22 = tau_33 = 800
void testFiniteStrainUniaxialResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  auto solver = makeSolver( matName, matProps );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 0.1; // F_11 = 1.1
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 1.0;
  step.dTMax                        = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  // For F = diag(1.1, 1, 1), Biot strain U-I = diag(0.1, 0, 0)
  // S_biot_11 = 24000*0.1 = 2400, S_biot_22 = S_biot_33 = 8000*0.1 = 800
  // PK2_11 = 2*2400/(1.1+1.1) = 2400/1.1
  // PK2_22 = PK2_33 = 2*800/2 = 800
  // tau_11 = 1.1^2 * 2400/1.1 = 1.1 * 2400 = 2640
  // tau_22 = tau_33 = 800
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 2640.0;
  stressTarget( 1, 1 ) = 800.0;
  stressTarget( 2, 2 ) = 800.0;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-8 ),
                           "I-3: Finite strain uniaxial response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-4: Stress tensor symmetry for an arbitrary deformation (isotropic)
void testStressTensorSymmetry()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  auto solver = makeSolver( matName, matProps );

  // Arbitrary non-symmetric displacement gradient
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 0.01;
  step.gradUIncrementTarget( 0, 1 ) = 0.06;
  step.gradUIncrementTarget( 0, 2 ) = -0.03;
  step.gradUIncrementTarget( 1, 0 ) = 0.06;
  step.gradUIncrementTarget( 1, 1 ) = 0.02;
  step.gradUIncrementTarget( 1, 2 ) = 0.04;
  step.gradUIncrementTarget( 2, 0 ) = -0.03;
  step.gradUIncrementTarget( 2, 1 ) = 0.04;
  step.gradUIncrementTarget( 2, 2 ) = -0.05;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 1.0;
  step.dTMax                        = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  throwExceptionOnFailure( checkIfEqual( finalStress, Fastor::transpose( finalStress ), 1e-10 ),
                           "I-4: Kirchhoff stress tensor symmetry check failed in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-5: Objectivity test - tau(Q*F) = Q * tau(F) * Q^T
// For a purely elastic isotropic material, rotating the spatial frame should
// transform the Kirchhoff stress accordingly.
void testObjectivity()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  // First compute stress for original deformation
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

  auto solver = makeSolver( matName, matProps );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = F - Marmot::FastorStandardTensors::Spatial3D::I;
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.timeStart                   = 0.0;
  step.timeEnd                     = 1.0;
  step.dTStart                     = 1.0;
  step.dTMax                       = 1.0;

  solver.addStep( step );
  solver.solve();

  Tensor33d tauUnrotated = solver.getHistory().back().stress;

  // Now test for several rotation angles about the z-axis
  for ( int phi_deg = 30; phi_deg <= 180; phi_deg += 30 ) {
    double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.0;

    // F_rotated = Q * F
    Tensor33d F_rotated = einsum< ik, kj, to_ij >( Q, F );

    auto solver_rot = makeSolver( matName, matProps );

    MarmotMaterialPointSolverFiniteStrain::Step step_rot;
    step_rot.gradUIncrementTarget        = F_rotated - Marmot::FastorStandardTensors::Spatial3D::I;
    step_rot.stressIncrementTarget       = Tensor33d( 0.0 );
    step_rot.isGradUComponentControlled  = Tensor33t< bool >( true );
    step_rot.isStressComponentControlled = Tensor33t< bool >( false );
    step_rot.timeStart                   = 0.0;
    step_rot.timeEnd                     = 1.0;
    step_rot.dTStart                     = 1.0;
    step_rot.dTMax                       = 1.0;

    solver_rot.addStep( step_rot );
    solver_rot.solve();

    Tensor33d tauNew = solver_rot.getHistory().back().stress;

    // Expected: tau(Q*F) = Q * tau(F) * Q^T
    Tensor33d tauExpected = einsum< iI, IJ, jJ, to_ij >( Q, tauUnrotated, Q );

    throwExceptionOnFailure( checkIfEqual( tauNew, tauExpected, 1e-8 ),
                             "I-5: Objectivity test failed for phi_deg=" + std::to_string( phi_deg ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-6: Orthotropic material - uniaxial stretch in direction 1
// Reference values computed from the Biot hyperelastic formulation for the orthotropic stiffness tensor.
// Material: E1=20000, E2=10000, E3=15000, nu12=0.25, nu13=0.35, nu23=0.30, G12=4000, G13=6000, G23=5000
// For F = diag(1.1, 1, 1): tau_11 ≈ 3555.30, tau_22 ≈ 1100.29, tau_33 ≈ 1461.32
void testOrthotropicUniaxialResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getOrthotropicProperties();

  auto solver = makeSolver( matName, matProps );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 0.1; // F_11 = 1.1
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 1.0;
  step.dTMax                        = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  // Reference values from analytical computation (Python verification)
  // C++ calls stiffnessTensor(E1,E2,E3, nu12, nu23, nu13, G12, G23, G13)
  // where nu23 and nu13 arguments correspond to mat[5]=0.30 and mat[4]=0.35 respectively
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 3.555300859598857e+03;
  stressTarget( 1, 1 ) = 1.100286532951290e+03;
  stressTarget( 2, 2 ) = 1.461318051575932e+03;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-6 ),
                           "I-6: Orthotropic uniaxial response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-7: Orthotropic material - simple shear deformation
// For F = I + 0.1 * e1 ⊗ e2 with the orthotropic material
void testOrthotropicSimpleShearResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getOrthotropicProperties();

  auto solver = makeSolver( matName, matProps );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 1 ) = 0.1; // F_12 = 0.1 (simple shear)
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 1.0;
  step.dTMax                        = 1.0;

  solver.addStep( step );
  solver.solve();

  auto history     = solver.getHistory();
  auto finalStress = history.back().stress;

  // Reference values from analytical computation (Python verification)
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 6.094122349057655e+01;
  stressTarget( 0, 1 ) = 4.009906851900923e+02;
  stressTarget( 1, 0 ) = 4.009906851900923e+02;
  stressTarget( 1, 1 ) = 2.015652085841339e+01;
  stressTarget( 2, 2 ) = 1.235906850394132e+01;

  throwExceptionOnFailure( checkIfEqual( finalStress, stressTarget, 1e-6 ),
                           "I-7: Orthotropic simple shear response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-8: Pure rigid body rotation should yield zero stress (isotropic)
// For F = Q (pure rotation), the Biot strain U - I = 0 since U = I for pure rotation.
void testPureRotationResponse()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicProperties();

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.0;

    auto solver = makeSolver( matName, matProps );

    MarmotMaterialPointSolverFiniteStrain::Step step;
    step.gradUIncrementTarget        = Q - Marmot::FastorStandardTensors::Spatial3D::I;
    step.stressIncrementTarget       = Tensor33d( 0.0 );
    step.isGradUComponentControlled  = Tensor33t< bool >( true );
    step.isStressComponentControlled = Tensor33t< bool >( false );
    step.timeStart                   = 0.0;
    step.timeEnd                     = 1.0;
    step.dTStart                     = 1.0;
    step.dTMax                       = 1.0;

    solver.addStep( step );
    solver.solve();

    Tensor33d finalStress = solver.getHistory().back().stress;

    throwExceptionOnFailure( checkIfEqual( finalStress, Tensor33d( 0.0 ), 1e-10 ),
                             "I-8: Pure rotation (phi_deg=" + std::to_string( phi_deg ) +
                               ") should yield zero stress in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-9: Viscoelastic stress relaxation
// Under a constant deformation, the Kirchhoff stress should relax from the instantaneous
// to the long-term value. With gamma=0.3, the long-term stress should be (1 - 0.3) = 0.7
// times the instantaneous elastic stress.
void testViscoelasticRelaxation()
{
  const std::string matName  = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto              matProps = getIsotropicViscoelasticProperties();

  auto solver = makeSolver( matName, matProps );

  // Step 1: Apply deformation instantaneously (very small time increment)
  MarmotMaterialPointSolverFiniteStrain::Step stepLoad;
  stepLoad.gradUIncrementTarget         = Tensor33d( 0.0 );
  stepLoad.gradUIncrementTarget( 0, 0 ) = 0.1; // F_11 = 1.1
  stepLoad.stressIncrementTarget        = Tensor33d( 0.0 );
  stepLoad.isGradUComponentControlled   = Tensor33t< bool >( true );
  stepLoad.isStressComponentControlled  = Tensor33t< bool >( false );
  stepLoad.timeStart                    = 0.0;
  stepLoad.timeEnd                      = 1e-4;
  stepLoad.dTStart                      = 1e-4;
  stepLoad.dTMax                        = 1e-4;

  solver.addStep( stepLoad );

  // Step 2: Hold deformation constant for a long time (>> tau=10, so almost fully relaxed)
  // After many relaxation times, stress ≈ (1 - gamma) * instantaneous
  MarmotMaterialPointSolverFiniteStrain::Step stepHold;
  stepHold.gradUIncrementTarget        = Tensor33d( 0.0 ); // no additional deformation
  stepHold.stressIncrementTarget       = Tensor33d( 0.0 );
  stepHold.isGradUComponentControlled  = Tensor33t< bool >( true );
  stepHold.isStressComponentControlled = Tensor33t< bool >( false );
  stepHold.timeStart                   = 1e-4;
  stepHold.timeEnd                     = 1000.0; // 100 * tau, so fully relaxed
  stepHold.dTStart                     = 10.0;
  stepHold.dTMax                       = 100.0;

  solver.addStep( stepHold );
  solver.solve();

  auto history       = solver.getHistory();
  auto stressRelaxed = history.back().stress;

  // Compute the instantaneous (elastic) stress for F = diag(1.1, 1, 1)
  // isotropic E=20000, nu=0.25: tau_11 = 2640, tau_22 = tau_33 = 800
  // After full relaxation: stress = (1 - gamma) * instantaneous = 0.7 * instantaneous
  const double gamma       = 0.3;
  const double longTermFac = 1.0 - gamma;

  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = longTermFac * 2640.0;
  stressTarget( 1, 1 ) = longTermFac * 800.0;
  stressTarget( 2, 2 ) = longTermFac * 800.0;

  throwExceptionOnFailure( checkIfEqual( stressRelaxed, stressTarget, 1.0 ),
                           "I-9: Viscoelastic relaxation failed - long-term stress does not match "
                           "expected (1-gamma)*instantaneous in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-10: Substepped variant produces same result as regular variant
// Both FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY and _SUBSTEPPED should give
// the same stress for a given deformation (within numerical tolerance).
void testSubsteppedVariantConsistency()
{
  auto matProps = getOrthotropicProperties();

  // Define a deformation increment
  Tensor33d gradUIncrement( 0.0 );
  gradUIncrement( 0, 0 ) = 0.05;
  gradUIncrement( 1, 1 ) = 0.02;
  gradUIncrement( 2, 2 ) = -0.01;

  auto createStep = [&]() {
    MarmotMaterialPointSolverFiniteStrain::Step step;
    step.gradUIncrementTarget        = gradUIncrement;
    step.stressIncrementTarget       = Tensor33d( 0.0 );
    step.isGradUComponentControlled  = Tensor33t< bool >( true );
    step.isStressComponentControlled = Tensor33t< bool >( false );
    step.timeStart                   = 0.0;
    step.timeEnd                     = 1.0;
    step.dTStart                     = 1.0;
    step.dTMax                       = 1.0;
    return step;
  };

  // Regular variant
  std::string matNameRegular = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY";
  auto        solverRegular  = makeSolver( matNameRegular, matProps );
  solverRegular.addStep( createStep() );
  solverRegular.solve();
  Tensor33d stressRegular = solverRegular.getHistory().back().stress;

  std::vector< double > matPropsSub;
  matPropsSub.push_back( 1.0 ); // one substep
  for ( auto& val : matProps ) {
    matPropsSub.push_back( val );
  }
  // Substepped variant
  std::string matNameSub = "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED";
  auto        solverSub  = makeSolver( matNameSub, matPropsSub );
  solverSub.addStep( createStep() );
  solverSub.solve();
  Tensor33d stressSub = solverSub.getHistory().back().stress;

  throwExceptionOnFailure( checkIfEqual( stressRegular, stressSub, 1e-8 ),
                           "I-10: Substepped variant gives different result than regular variant in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,             // I-1: F=I gives zero stress
    testSmallStrainUniaxialResponse,    // I-2: Small strain uniaxial (isotropic)
    testFiniteStrainUniaxialResponse,   // I-3: Finite strain uniaxial (isotropic)
    testStressTensorSymmetry,           // I-4: Stress symmetry for arbitrary deformation
    testObjectivity,                    // I-5: Objectivity (tau(Q*F) = Q*tau(F)*Q^T)
    testOrthotropicUniaxialResponse,    // I-6: Orthotropic uniaxial stretch
    testOrthotropicSimpleShearResponse, // I-7: Orthotropic simple shear
    testPureRotationResponse,           // I-8: Pure rotation gives zero stress
    testViscoelasticRelaxation,         // I-9: Viscoelastic stress relaxation
    testSubsteppedVariantConsistency,   // I-10: Substepped == regular
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
