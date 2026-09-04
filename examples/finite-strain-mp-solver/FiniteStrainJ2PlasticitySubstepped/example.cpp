#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
// Ensure the registration unit is linked or included if using static registration
// #include "Marmot/MaterialRegistration.cpp"

int main()
{
  using namespace Marmot::Solvers;
  using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;

  // 1) Define the material model
  std::string materialName = "FINITESTRAINJ2PLASTICITY_SUBSTEPPED";

  // Define Properties
  double properties[] = {
    20,   // nSubsteps (20 substeps per global increment)
    35,   // K (Bulk Modulus)
    15,   // G (Shear Modulus)
    0.45, // fy (Initial Yield Stress)
    0.75, // fyInf (Saturated Yield Stress)
    0.5,  // eta (Saturation rate)
    0.1,  // H (Linear Hardening)
    1.0,  // implementationType (1 = Full Analytical Return Mapping)
    0.0   // density
  };
  int nProps = 9;

  // 2) Configure solver options
  MarmotMaterialPointSolverFiniteStrain::SolverOptions options;
  options.maxIterations       = 15;
  options.residualTolerance   = 1e-8;
  options.correctionTolerance = 1e-8;

  // 3) Create solver instance
  MarmotMaterialPointSolverFiniteStrain solver( materialName, properties, nProps, options );

  // 5) Define Loading Steps

  // --- Step 1: Loading (Uniaxial extension to 100%) ---
  MarmotMaterialPointSolverFiniteStrain::Step step1;
  step1.timeStart = 0.0;
  step1.timeEnd   = 1.0;
  step1.dTStart   = 1e-0;
  step1.dTMin     = 1e-1;
  step1.dTMax     = 0.5;

  step1.gradUIncrementTarget        = Tensor33d( 0.0 );
  step1.stressIncrementTarget       = Tensor33d( 0.0 );
  step1.isGradUComponentControlled  = Tensor33t< bool >( false );
  step1.isStressComponentControlled = Tensor33t< bool >( true );

  // Control gradU_11 (stretch to 10%)
  step1.gradUIncrementTarget( 0, 0 )        = 1;
  step1.isGradUComponentControlled( 0, 0 )  = true;
  step1.isStressComponentControlled( 0, 0 ) = false;

  solver.addStep( step1 );

  // --- Step 2: Unloading (Return to 0% global strain) ---
  MarmotMaterialPointSolverFiniteStrain::Step step2 = step1;
  step2.timeStart                                   = 1.0;
  step2.timeEnd                                     = 2.0;

  step2.gradUIncrementTarget( 0, 0 ) = -1;

  solver.addStep( step2 );

  // 6) Run the solver
  try {
    solver.solve();
  }
  catch ( const std::exception& e ) {
    std::cerr << "\n!!! SOLVER FAILED !!!\nError: " << e.what() << std::endl;
    return 1;
  }

  // 7) Output the results
  solver.exportHistoryToCSV( "j2_substepped_results.csv" );

  return 0;
}
