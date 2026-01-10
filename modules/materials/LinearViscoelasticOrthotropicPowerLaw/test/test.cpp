#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot::Testing;
using namespace Marmot::Solvers;

Eigen::Vector< double, 23 > getMaterialPropertiesIsotropic()
{
  double E  = 2e5;
  double nu = 0.2;
  double G  = E / ( 2 * ( 1 + nu ) );
  // material properties
  std::vector< double > materialProperties;
  // Young's moduli
  materialProperties.push_back( 1.0 );
  materialProperties.push_back( E );
  materialProperties.push_back( E );
  materialProperties.push_back( E );
  // Poisson's ratios
  materialProperties.push_back( nu );
  materialProperties.push_back( nu );
  materialProperties.push_back( nu );
  // shear moduli
  materialProperties.push_back( G );
  materialProperties.push_back( G );
  materialProperties.push_back( G );
  // viscoelastic parameters
  materialProperties.push_back( 0.5 );
  materialProperties.push_back( 0.1 );
  materialProperties.push_back( 2 );
  materialProperties.push_back( 10 );
  materialProperties.push_back( 0.0001 );
  materialProperties.push_back( 3.1622776601683795 );
  materialProperties.push_back( 1.0 );
  // coordinate system
  materialProperties.push_back( 1.0 );
  materialProperties.push_back( 0.5 );
  materialProperties.push_back( 0.0 );
  materialProperties.push_back( -0.5 );
  materialProperties.push_back( 1.0 );
  materialProperties.push_back( 1.0 );

  return Eigen::Vector< double, 23 >( &materialProperties[0] );
}

void testLinearViscoelasticOrthotropicPowerLawIsotropic()
{

  auto        materialProperties = getMaterialPropertiesIsotropic();
  auto        solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string matName            = "LINEARVISCOELASTICORTHOTROPICPOWERLAW";
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

int main()
{
  testLinearViscoelasticOrthotropicPowerLawIsotropic();
  return 0;
}
