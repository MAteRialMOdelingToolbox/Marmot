#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"

int main()
{

  double appliedStrain = 0.1;

  using namespace Marmot::Materials;
  using namespace Marmot::Math::Tensors::FastorTensors;

  // 1) Define the material model
  // std::string materialName = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED";
  std::string materialName = "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY";
  // layout: [hyperelasticBase, onlyShearCreep, <elastic params of the base>, nMaxwell, (gamma_i, tau_i)...]
  // here: NeoHooke (0) with K=3500, G=1500, full creep (0), two Maxwell elements
  double properties[] = { 0, 0, 3500, 1500, 2, 0.5, 10000000, 0.25, 10 };
  int    nProps       = 9;

  // 2) Configure solver options
  MaterialPointSolverFiniteStrain::SolverOptions options;
  options.maxIterations       = 25;
  options.residualTolerance   = 1e-10;
  options.correctionTolerance = 1e-10;

  // 3) Create solver instance
  MaterialPointSolverFiniteStrain solver( materialName, properties, nProps, options );

  // 5) Define a loading step (Uniaxial extension in 11-direction)
  MaterialPointSolverFiniteStrain::Step step;
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
  MaterialPointSolverFiniteStrain::Step relaxStep;
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
