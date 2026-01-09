#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"

int main()
{
  using namespace Marmot;
  using namespace Fastor;

  // 1) Define the material model
  std::string materialName = "COMPRESSIBLENEOHOOKE";
  double      properties[] = { 3500, 1500 }; // Example: Bulk and Shear moduli
  int         nProps       = 2;

  // 2) Configure solver options
  MarmotMaterialPointSolverFiniteStrain::SolverOptions options;
  options.maxIterations       = 25;
  options.residualTolerance   = 1e-10;
  options.correctionTolerance = 1e-10;

  // 3) Create solver instance
  MarmotMaterialPointSolverFiniteStrain solver( materialName, properties, nProps, options );

  // 4) Set the initial state
  int nSV = 0;
  solver.getNumberOfStateVariables( nSV );
  Eigen::VectorXd initialSV = Eigen::VectorXd::Zero( nSV );

  // Initial stress is zero
  solver.setInitialState( Tensor33d( 0.0 ), initialSV );

  // 5) Define a loading step (Uniaxial extension in 11-direction)
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.timeStart = 0.0;
  step.timeEnd   = 1.0;
  step.dTStart   = 0.1;
  step.dTMin     = 1e-6;
  step.dTMax     = 0.5;

  // Initialize targets and flags
  step.gradUIncrementTarget        = Tensor33d( 0.0 );
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( false );
  step.isStressComponentControlled = Tensor33t< bool >( true );

  // Setup Control:
  // - Control gradU_11 (stretch to 50%)
  step.gradUIncrementTarget( 0, 0 )       = 0.5;
  step.isGradUComponentControlled( 0, 0 ) = true;

  // - Unset control on tau_11
  step.isStressComponentControlled( 0, 0 ) = false;

  // 6) Run the solver
  solver.addStep( step );
  solver.solve();

  // 7) Output the results
  solver.exportHistoryToCSV( "finite_strain_history.csv" );

  return 0;
}
