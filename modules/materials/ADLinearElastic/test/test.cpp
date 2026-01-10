#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot::Testing;
using namespace Marmot::Solvers;

auto initializeSolverADLinearElastic()
{
  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  return std::unique_ptr< MarmotMaterialPointSolverHypoElastic >(
    new MarmotMaterialPointSolverHypoElastic( matName,
                                              &materialProperties[0],
                                              materialProperties.size(),
                                              solveropts ) );
}

void testADLinearElasticObjectivity()
{

  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;

  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };

  step.dTStart = 1.0;

  solver.addStep( step );
  // check coordinate invariance
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << "Turbokreisel failed!" );
}

void testADLinearElastic()
{

  std::vector< double > materialProperties = { 210000., 0.3 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "ADLINEARELASTIC";
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;

  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };

  step.dTStart = 1.0;

  solver.addStep( step );

  // solve
  solver.solve();

  // reference solution
  const Marmot::Vector6d stressTarget =
    { 1627.90061538, 416.523692308, -864.896307692, 11.0128846154, 129.673384615, 462.666346154 };

  auto history = solver.getHistory();

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << "comparison with reference solution failed" );
}

int main()
{

  std::vector< std::function< void( void ) > > tests = {
    testADLinearElastic,
    testADLinearElasticObjectivity,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
