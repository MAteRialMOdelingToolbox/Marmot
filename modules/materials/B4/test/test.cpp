#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot::Testing;

Eigen::Vector< double, 22 > getMaterialPropertiesB4()
{
  // material properties
  Eigen::Vector< double, 22 > materialProperties;
  // elastic parameters
  double nu = 0.2;
  double q1 = 20.6;
  // basic creep parameters
  double q2           = 91;
  double q3           = 4.80;
  double q4           = 5.9;
  double n            = 0.1;
  double m            = 0.5;
  int    nKelvinBasic = 12;
  double minTauBasic  = 1e-5;
  // autogenous shrinkage parameters
  double ultimateAutogenousShrinkageStrain = 0.; // -0.0001
  double autogenousShrinkageHalfTime       = 3.;
  double alpha                             = 1.45;
  double rt                                = -4.5;
  // drying shrinkage parameters
  double ultimateDryingShrinkageStrain = -0.0015;
  double dryingShrinkageHalfTime       = 90.;
  double dryingStart                   = 7.;
  double hEnv                          = 1.;
  // drying creep paramters
  double q5            = 400;
  int    nKelvinDrying = 11;
  double minTauDrying  = 2e-4;
  // additional parameters
  double castTime   = -100;
  double timeToDays = 1;
  materialProperties << nu, q1, q2, q3, q4, n, m, nKelvinBasic, minTauBasic, ultimateAutogenousShrinkageStrain,
    autogenousShrinkageHalfTime, alpha, rt, ultimateDryingShrinkageStrain, dryingShrinkageHalfTime, dryingStart, hEnv,
    q5, nKelvinDrying, minTauDrying, castTime, timeToDays;

  return materialProperties;
}

void testB4()
{
  auto        materialProperties = getMaterialPropertiesB4();
  auto        solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName            = "B4";
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

  // expected stress
  Marmot::Vector6d stressTarget;
  stressTarget << 0.377092, 0.377092, 0.150837, 0., 0., 0.;

  // read history
  auto history = solver.getHistory();
  // get computed stress
  Marmot::Vector6d stress = history.back().stress;

  std::cout << "Computed stress: " << stress.transpose() << std::endl;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  Marmot::Testing::throwExceptionOnFailure( Marmot::Testing::checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                                            "Stress computation for " + std::string( __PRETTY_FUNCTION__ ) );
}

void testB4CoordinateInvariance()
{

  auto        materialProperties = getMaterialPropertiesB4();
  auto        solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName            = "B4";
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
    testB4,
    // testB4CoordinateInvariance
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
