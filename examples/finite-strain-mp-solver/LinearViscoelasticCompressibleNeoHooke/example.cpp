#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"

int main( int argc, char** argv )
{

  // get applied strain from command line
  // double appliedStrain = argv[1] ? atof( argv[1] ) : 0.5;
  // int    nMaxwell      = argv[2] ? atoi( argv[2] ) : 2;

  double appliedStrain = 0.1;

  using namespace Marmot::Solvers;
  using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;

  // 1) Define the material model
  // std::string materialName = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED";
  std::string materialName = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  // double      properties[] = { 2, 750, 75, 2./3500,  double( nMaxwell ), 0.5, 1, 0.25, 10 };
  // double properties[] = { 1, 0, 3500,1500, 1, 0.5, 1, 0.25, 10 };
  double properties[] = { 0, 3500, 1500, 2, 0.5, 10000000, 0.25, 10 };
  // double      properties[] = { 0, 3500, 1500,  double( nMaxwell ), 0.5, 1, 0.25, 10 };
  int nProps = 8;

  // 2) Configure solver options
  MarmotMaterialPointSolverFiniteStrain::SolverOptions options;
  options.maxIterations       = 25;
  options.residualTolerance   = 1e-10;
  options.correctionTolerance = 1e-10;

  // 3) Create solver instance
  MarmotMaterialPointSolverFiniteStrain solver( materialName, properties, nProps, options );

  // 5) Define a loading step (Uniaxial extension in 11-direction)
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.timeStart = 0.0;
  step.timeEnd   = 1e-4;
  step.dTStart   = 1e-5;
  step.dTMin     = 1e-5;
  step.dTMax     = 1e-5;

  // Initialize targets and flags
  step.gradUIncrementTarget        = Tensor33d( 0.0 );
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( false );
  step.isStressComponentControlled = Tensor33t< bool >( true );

  // Setup Control:
  // - Control gradU_11 (stretch to 50%)
  step.gradUIncrementTarget( 0, 0 )       = appliedStrain;
  step.isGradUComponentControlled( 0, 0 ) = true;

  // - Unset control on tau_11
  step.isStressComponentControlled( 0, 0 ) = false;

  solver.addStep( step );

  // add a relaxation step
  MarmotMaterialPointSolverFiniteStrain::Step relaxStep;
  relaxStep.timeStart                          = 1e-9;
  relaxStep.timeEnd                            = 20;
  relaxStep.dTStart                            = .1;
  relaxStep.dTMin                              = 1e-6;
  relaxStep.dTMax                              = 0.5;
  relaxStep.maxIncrements                      = 1000;
  relaxStep.gradUIncrementTarget               = Tensor33d( 0.0 );
  relaxStep.stressIncrementTarget              = Tensor33d( 0.0 );
  relaxStep.isGradUComponentControlled         = Tensor33t< bool >( false );
  relaxStep.isStressComponentControlled        = Tensor33t< bool >( true );
  relaxStep.gradUIncrementTarget( 0, 0 )       = 0.0;
  relaxStep.isGradUComponentControlled( 0, 0 ) = true;

  // - Unset control on tau_11
  relaxStep.isStressComponentControlled( 0, 0 ) = false;

  // solver.addStep( relaxStep );

  solver.solve();

  // 7) Output the results
  solver.exportHistoryToCSV( "finite_strain_history.csv" );

  return 0;
}
