#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot::Testing;

Eigen::Vector< double, 7 > getMaterialPropertiesLinearViscoelasticPowerLaw()
{
  // material properties
  Eigen::Vector< double, 7 > materialProperties;
  // elastic parameters
  double E  = 2e5;
  double nu = 0.2;
  // viscoelastic parameters
  double m          = 0.5;
  double n          = 0.1;
  int    nKelvin    = 10;
  double minTau     = 0.0001;
  double timeToDays = 1.;
  materialProperties << E, nu, m, n, nKelvin, minTau, timeToDays;

  return materialProperties;
}

void testLinearViscoelasticPowerLaw()
{

  auto        materialProperties = getMaterialPropertiesLinearViscoelasticPowerLaw();
  auto        solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName            = "LINEARVISCOELASTICPOWERLAW";
  auto        solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // step 1: advance time to 28 days
  MarmotMaterialPointSolverHypoElastic::Step step1;

  step1.isStrainComponentControlled = { true, true, true, true, true, true };
  step1.isStressComponentControlled = { false, false, false, false, false, false };
  step1.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step1.strainIncrementTarget       = { 0.0, 0., 0., 0., 0.0, 0.0 };
  step1.timeEnd                     = 28.0;
  step1.dTStart                     = 28;

  solver.addStep( step1 );

  // step 1: apply strain
  MarmotMaterialPointSolverHypoElastic::Step step2;
  step2.isStrainComponentControlled = { true, true, true, true, true, true };
  step2.isStressComponentControlled = { false, false, false, false, false, false };
  step2.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step2.strainIncrementTarget       = { 10e-6, 10e-6, 0., 0., 0., 0. };
  step2.dTStart                     = 0.01;
  step2.timeStart                   = 28.0;
  step2.timeEnd                     = 28.01;

  solver.addStep( step2 );

  // step 2: hold strain constant
  MarmotMaterialPointSolverHypoElastic::Step step3;
  step3.isStrainComponentControlled = { true, true, true, true, true, true };
  step3.isStressComponentControlled = { false, false, false, false, false, false };
  step3.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step3.strainIncrementTarget       = { 0., 0., 0., 0., 0., 0. };
  step3.dTStart                     = 100;
  step3.timeStart                   = 28.01;
  step3.timeEnd                     = 128.01;

  solver.addStep( step3 );

  solver.solve();

  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  Marmot::Vector6d stressTarget;
  stressTarget << 2.41296e-05, 2.41296e-05, 9.65182e-06, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  Marmot::Testing::throwExceptionOnFailure( Marmot::Testing::checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                                            "Stress computation failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

void testLinearViscoelasticPowerLawCoordinateInvariance()
{

  auto        materialProperties = getMaterialPropertiesLinearViscoelasticPowerLaw();
  auto        solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName            = "LINEARVISCOELASTICPOWERLAW";
  auto        solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;

  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0., 0., 0., 0.03, 0.0 };

  step.dTStart = 1.0;

  solver.addStep( step );

  // check coordinate invariance
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

int main()
{
  std::vector< std::function< void() > > tests = {
    testLinearViscoelasticPowerLaw,
    testLinearViscoelasticPowerLawCoordinateInvariance,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
