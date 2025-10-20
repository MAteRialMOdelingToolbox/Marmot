#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"

int main()
{

  double                                              matProps[6] = { 210e3, 0.3, 200, 2100, 20, 20 };
  const int                                           nProps      = 6;
  std::string                                         name        = "ADVONMISES";
  MarmotMaterialPointSolverHypoElastic::SolverOptions options;

  auto solver = MarmotMaterialPointSolverHypoElastic( name, &matProps[0], nProps, options );

  MarmotMaterialPointSolverHypoElastic::Step step;

  step.stressIncrementTarget       = { 300.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.isStressComponentControlled = { true, true, true, false, false, false };
  step.isStrainComponentControlled = { false, false, false, true, true, true };
  step.dTStart                     = 0.01;

  solver.addStep( step );
  solver.addStep( step );

  solver.solve();

  auto history = solver.getHistory();

  solver.printHistory();

  solver.exportHistoryToCSV( "MarmotMaterialPointSolverResults.csv" );

  return 0;
}
