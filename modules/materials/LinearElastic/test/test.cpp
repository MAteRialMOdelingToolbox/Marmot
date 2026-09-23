#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>

// Use namespaces for brevity
using namespace Marmot::Testing;
using namespace Marmot::Solvers;

// Function to test the isotropic material response for a normal strain increment
void testMaterialResponse()
{
  // Define material parameters (Young's modulus and Poisson's ratio)
  std::vector< double > materialProperties = { 20000, 0.25 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.001, 0., 0., 0., 0.0, 0.0 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget;
  stressTarget << 24., 8., 8., 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                           "Stress computation failed for normal strain in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the isotropic material response for a shear strain increment
void testShearMaterialResponse()
{

  // Define material parameters (Young's modulus and Poisson's ratio)
  std::vector< double > materialProperties = { 20000, 0.25 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.0, 0., 0., 0.001, 0.001, 0.001 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear;
  stressTargetShear << 0., 0., 0., 8., 8., 8.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear, 1e-10 ),
                           "Stress computation failed for shear strain in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the transversely isotropic material response for a normal strain increment
void testTransverseIsotropicMaterialResponse()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  std::vector< double > materialProperties = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.001, 0., 0., 0., 0.0, 0.0 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_transverseIsotropic;
  stressTarget_transverseIsotropic << 21.96078431372549, 3.9215686274509802, 3.9215686274509802, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_transverseIsotropic, 1e-10 ),
                           "Stress computation failed for transverse isotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the transversely isotropic material response for a shear strain increment
void testTransverseIsotropicShearMaterialResponse()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  std::vector< double > materialProperties = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 1, 0 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.0, 0., 0., 0.001, 0.001, 0.001 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  // Define the expected stress values for the applied shear strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear_transverseIsotropic;
  stressTargetShear_transverseIsotropic << 0., 0., 0., 4., 4., 3.846153846153847;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear_transverseIsotropic, 1e-10 ),
                           "Stress computation failed for transverse isotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a normal strain increment
void testOrthotropicMaterialResponse()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  std::vector< double > materialProperties =
    { 20000., 10000., 15000., 0.25, 0.3, 0.35, 4000, 5000, 6000, 1, 0, 0, 0, 1, 0 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.001, 0., 0., 0., 0.0, 0.0 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;
  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_orthotropic;
  stressTarget_orthotropic << 24.62633451957296, 5.800711743772242, 9.074733096085410, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_orthotropic, 1e-10 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a shear strain increment
void testOrthotropicShearMaterialResponse()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  std::vector< double > materialProperties =
    { 20000, 10000, 15000, 0.25, 0.3, 0.35, 4000, 5000, 6000, 1, 0, 0, 0, 1, 0 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.0, 0., 0., 0.001, 0.001, 0.001 };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;
  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear_orthotropic;
  stressTargetShear_orthotropic << 0., 0., 0., 4., 6., 5.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear_orthotropic, 1e-10 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a normal strain increment
void testOrthotropicMaterialResponseRotation()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  std::vector< double > materialProperties =
    { 1000, 30, 30, 0.009, 0, 0, 50, 50, 50, 0.906307787, 0.422618262, 0, -0.422618262, 0.906307787, 0 };

  auto        solveropts = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName    = "LINEARELASTIC";
  auto        solver     = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define step: apply normal strain increment
  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { -2.4e-4, -3.9e-4, 0., 1.76e-3, 0., 0. };

  // add step to solver
  solver.addStep( step );
  // solve
  solver.solve();
  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  // Define the expected st ress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_orthotropic;
  stressTarget_orthotropic << 2.898786736584660e-01, 8.616047561468707e-02, 0., 2.004527953677034e-01, 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_orthotropic, 1e-8 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for isotropic material
void testGetDensityIsotropic()
{
  // Define material parameters (Young's modulus, Poisson's ratio and Density)
  const double materialProperties[3] = { 20000, 0.25, 10 };
  const int    nMaterialProperties   = 3;

  auto mat = Marmot::Materials::LinearElastic( materialProperties, nMaterialProperties, 1 );

  // Defined density at start (materialProperties[3])
  double expectedDensityI = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat.getDensity( nullptr ), expectedDensityI, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for transversely isotropic material
void testGetDensityTransverselyIsotropic()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  const double materialProperties[12] = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1, 10 };
  const int    nMaterialProperties    = 12;

  auto mat = Marmot::Materials::LinearElastic( materialProperties, nMaterialProperties, 1 );

  // Defined density at start (materialProperties[3])
  double expectedDensityTI = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat.getDensity( nullptr ), expectedDensityTI, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for transversely isotropic material
void testGetDensityOrthotropic()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  const double materialProperties[16] =
    { 1000, 30, 30, 0.009, 0, 0, 50, 0, 0, 0.906307787, 0.422618262, 0, -0.422618262, 0.906307787, 0, 10 };
  const int nMaterialProperties = 16;

  auto mat = Marmot::Materials::LinearElastic( materialProperties, nMaterialProperties, 1 );

  // Defined density at start (materialProperties[3])
  double expectedDensityO = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat.getDensity( nullptr ), expectedDensityO, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

// ─────────────────────────────────────────────────────────────────────────────
// The following tests exercise MarmotMaterialHypoElastic's default (base-class)
// implementations, none of which LinearElastic overrides: computeStressExplicit,
// computePlaneStress, computeUniaxialStress, getMaximumWaveSpeed, getStateView,
// initializeYourself and setCharacteristicElementLength.
// ─────────────────────────────────────────────────────────────────────────────

namespace {
  std::array< double, 2 >          isotropicProps_ = { 20000., 0.25 };
  Marmot::Materials::LinearElastic makeIsotropicMaterial()
  {
    return { isotropicProps_.data(), 2, 1 };
  }
} // namespace

void testComputeStressExplicitMatchesComputeStress()
{
  using namespace Marmot;
  auto matA = makeIsotropicMaterial();
  auto matB = makeIsotropicMaterial();

  Vector6d dStrain = Vector6d::Zero();
  dStrain( 0 )     = 1e-3;
  dStrain( 3 )     = 2e-4;

  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0., 1. };

  MarmotMaterialHypoElastic::state3D stateA;
  Matrix6d                           tanA;
  // LinearElastic re-declares computeStress() under `protected`, so call it through the base
  // class (as production code always does, via a MarmotMaterialHypoElastic*).
  MarmotMaterialHypoElastic& matABase = matA;
  matABase.computeStress( stateA, tanA, dStrain, timeInfo );

  MarmotMaterialHypoElastic::state3D stateB;
  matB.computeStressExplicit( stateB, dStrain, timeInfo );

  throwExceptionOnFailure( checkIfEqual< double >( stateA.stress, stateB.stress, 1e-12 ),
                           "computeStressExplicit() must match computeStress() ignoring the tangent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputePlaneStressMatchesClassicalPlaneStressLaw()
{
  using namespace Marmot;
  auto mat = makeIsotropicMaterial();

  const double E = isotropicProps_[0], nu = isotropicProps_[1];

  Vector3d dStrain2D;
  dStrain2D << 1e-3, 0.5e-3, 2e-4; // eps_xx, eps_yy, gamma_xy

  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0., 1. };
  MarmotMaterialHypoElastic::state2D  state2D;
  Matrix3d                            tan2D;

  mat.computePlaneStress( state2D, tan2D, dStrain2D, timeInfo );

  // Classical isotropic plane-stress law (sigma_zz = 0, both eps_xx and eps_yy prescribed):
  //   sigma_xx = E/(1-nu^2) * (eps_xx + nu*eps_yy)
  //   sigma_yy = E/(1-nu^2) * (eps_yy + nu*eps_xx)
  //   sigma_xy = G * gamma_xy
  const double planeStressModulus = E / ( 1. - nu * nu );
  const double G                  = E / ( 2. * ( 1. + nu ) );
  Vector3d     expected;
  expected << planeStressModulus * ( dStrain2D( 0 ) + nu * dStrain2D( 1 ) ),
    planeStressModulus * ( dStrain2D( 1 ) + nu * dStrain2D( 0 ) ), G * dStrain2D( 2 );

  throwExceptionOnFailure( checkIfEqual< double >( state2D.stress, expected, 1e-6 ),
                           "computePlaneStress() does not match the classical isotropic plane-stress law in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputeUniaxialStressMatchesE()
{
  using namespace Marmot;
  auto mat = makeIsotropicMaterial();

  const double E         = isotropicProps_[0];
  const double dStrain1D = 1e-3;

  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0., 1. };
  MarmotMaterialHypoElastic::state1D  state1D;
  double                              tan1D;

  mat.computeUniaxialStress( state1D, tan1D, dStrain1D, timeInfo );

  throwExceptionOnFailure( checkIfEqual( state1D.stress, E * dStrain1D, 1e-6 ),
                           "computeUniaxialStress() does not reduce to sigma = E*eps in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( tan1D, E, 1e-6 ),
                           "computeUniaxialStress() tangent does not reduce to E in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testGetMaximumWaveSpeedMatchesPWaveModulus()
{
  using namespace Marmot;
  const std::array< double, 3 > props = { 20000., 0.25, 2400. };
  auto                          mat   = Marmot::Materials::LinearElastic( props.data(), 3, 1 );

  const double E = props[0], nu = props[1], rho = props[2];
  const double mu       = E / ( 2. * ( 1. + nu ) );
  const double lambda   = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  const double expected = std::sqrt( ( lambda + 2. * mu ) / rho );

  MarmotMaterialHypoElastic::state3D state;
  const double                       waveSpeed = mat.getMaximumWaveSpeed( state );

  throwExceptionOnFailure( checkIfEqual( waveSpeed, expected, 1e-5 ),
                           "getMaximumWaveSpeed() does not match sqrt((lambda+2*mu)/rho) in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testGetStateViewThrowsForMaterialWithNoStateVars()
{
  auto mat = makeIsotropicMaterial();

  bool threw = false;
  try {
    mat.getStateView( "anything", nullptr );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "getStateView() must throw for a material with no registered state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testInitializeYourselfZeroesStateVars()
{
  auto                  mat = makeIsotropicMaterial();
  std::vector< double > stateVars( 3, 42.0 );
  mat.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

  throwExceptionOnFailure( std::all_of( stateVars.begin(), stateVars.end(), []( double v ) { return v == 0.0; } ),
                           "initializeYourself() must zero all state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSetCharacteristicElementLengthIsStored()
{
  auto mat = makeIsotropicMaterial();
  mat.setCharacteristicElementLength( 0.5 );
  throwExceptionOnFailure( mat.characteristicElementLength == 0.5,
                           "setCharacteristicElementLength() must store the given length in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testMaterialResponse,                         // test for normal strain
    testShearMaterialResponse,                    // test for shear strain
    testTransverseIsotropicMaterialResponse,      // test for transverse isotropic normal strain
    testTransverseIsotropicShearMaterialResponse, // test for transverse isotropic shear strain
    testOrthotropicMaterialResponse,              // test for orthotropic normal strain
    testOrthotropicShearMaterialResponse,         // test for orthotropic shear strain
    testOrthotropicMaterialResponseRotation,      // test for orthotropic normal strain with rotation
    testGetDensityIsotropic,                      // test for density retrieval for isotropic case
    testGetDensityTransverselyIsotropic,          // test for density retrieval for transversely isotropic case
    testGetDensityOrthotropic,                    // test for density retrieval for orthotropic case
    testComputeStressExplicitMatchesComputeStress,
    testComputePlaneStressMatchesClassicalPlaneStressLaw,
    testComputeUniaxialStressMatchesE,
    testGetMaximumWaveSpeedMatchesPWaveModulus,
    testGetStateViewThrowsForMaterialWithNoStateVars,
    testInitializeYourselfZeroesStateVars,
    testSetCharacteristicElementLengthIsStored,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
