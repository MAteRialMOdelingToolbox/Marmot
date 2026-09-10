#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotTesting.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Solvers;
using namespace Marmot::FastorStandardTensors;

namespace {

  // A material that always fails to update its stress, used to exercise solveStep()'s
  // cutback-retry-then-give-up logic deterministically (no real material fails unconditionally).
  class AlwaysFailingFiniteStrainMaterial : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& ) const override
    {
      throw Marmot::StressUpdateFailed( "test-only material always fails" );
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  // A material whose Kirchhoff stress output is constant (independent of deformation) with a
  // singular (zero) tangent, used to force a NaN out of the Newton-loop's linear solve.
  class SingularTangentMaterial : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&,
                        const TimeIncrement& ) const override
    {
      response.tau                  = Tensor33d( 0.0 );
      response.tau( 0, 0 )          = 5.0;
      response.elasticEnergyDensity = 0.0;
      response.dissipation          = 0.0;
      tangents.dTau_dF              = Tensor3333d( 0.0 ); // singular
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  // A material whose Kirchhoff stress output is constant (independent of deformation), but with a
  // well-conditioned (identity) tangent, so the Newton loop computes well-defined, non-NaN
  // corrections every iteration without the residual ever actually shrinking.
  class WellConditionedNeverConvergingMaterial : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&,
                        const TimeIncrement& ) const override
    {
      response.tau                  = Tensor33d( 0.0 );
      response.tau( 0, 0 )          = 5.0;
      response.elasticEnergyDensity = 0.0;
      response.dissipation          = 0.0;

      tangents.dTau_dF = Tensor3333d( 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          tangents.dTau_dF( i, j, i, j ) = 1.0; // identity operator: well-conditioned, invertible
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  // A material that fails only for a genuinely nonzero time increment. solveStep() always
  // resolves its very first increment with dT == 0 (a trivial pass before dTStart is assigned on
  // the first successful increment), so this succeeds once before it starts failing -- letting a
  // subsequent retry actually halve an already-nonzero dT, rather than immediately hitting the
  // "dT already at or below dTMin" exhaustion check with the still-zero initial dT.
  class FailsForNonzeroTimeIncrementMaterial : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&,
                        const TimeIncrement& timeIncrement ) const override
    {
      if ( timeIncrement.dT > 0.0 )
        throw Marmot::StressUpdateFailed( "test-only material fails for any nonzero time increment" );
      response.tau                  = Tensor33d( 0.0 );
      response.elasticEnergyDensity = 0.0;
      response.dissipation          = 0.0;
      tangents.dTau_dF              = Tensor3333d( 0.0 );
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  const bool alwaysFailingRegistered = MarmotLibrary::MarmotMaterialFiniteStrainFactory::registerMaterial<
    AlwaysFailingFiniteStrainMaterial >( "TESTALWAYSFAILINGFINITESTRAIN" );
  const bool failsForNonzeroDTRegistered = MarmotLibrary::MarmotMaterialFiniteStrainFactory::registerMaterial<
    FailsForNonzeroTimeIncrementMaterial >( "TESTFAILSFORNONZERODTFINITESTRAIN" );
  const bool singularTangentRegistered = MarmotLibrary::MarmotMaterialFiniteStrainFactory::registerMaterial<
    SingularTangentMaterial >( "TESTSINGULARTANGENTFINITESTRAIN" );
  const bool neverConvergingRegistered = MarmotLibrary::MarmotMaterialFiniteStrainFactory::registerMaterial<
    WellConditionedNeverConvergingMaterial >( "TESTNEVERCONVERGINGFINITESTRAINSOLVER" );

} // namespace

// ---------------------------------------------------------------------------------------------
// setInitialState() / resetToInitialState()
// ---------------------------------------------------------------------------------------------
void testResetToInitialStateUndoesAccumulatedDeformation()
{
  // CompressibleNeoHooke is hyperelastic (path-independent): its stress is a pure function of the
  // current deformation gradient, not of the solver's "previous stress" -- so resetToInitialState()
  // is verified here by its actual load-bearing effect, undoing the internally accumulated
  // displacement gradient (gradU), rather than by an (inapplicable) stress-persistence check.
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-3;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );

  solver.addStep( step );
  solver.solve();
  const Tensor33d stressAfterFirstSolve = solver.getHistory().back().stress;

  solver.resetToInitialState();
  throwExceptionOnFailure( solver.getHistory().empty(),
                           "resetToInitialState() must clear the recorded history in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // resetToInitialState() does not clear queued steps, so avoid re-solving the (already-queued)
  // first step on top of itself.
  solver.clearSteps();
  solver.addStep( step );
  solver.solve();
  const Tensor33d stressAfterSecondSolve = solver.getHistory().back().stress;

  throwExceptionOnFailure( checkIfEqual( stressAfterFirstSolve( 0, 0 ), stressAfterSecondSolve( 0, 0 ), 1e-10 ),
                           "resetToInitialState() did not undo the internally accumulated displacement gradient "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSetInitialStateSeedsStateVariables()
{
  // Unlike stress (recomputed fresh from the deformation gradient for a hyperelastic law), state
  // variables genuinely persist: use a plastic material and verify that a custom initial hardening
  // variable is what ends up in the recorded history after a (near-)elastic step that doesn't
  // itself change it.
  std::vector< double > materialProperties = { 175000., 80800., 260., 580., 9., 70., 3. };
  std::string           matName            = "FINITESTRAINJ2PLASTICITY";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  const int nStateVars = solver.getNumberOfStateVariables();
  throwExceptionOnFailure( nStateVars == 10, // Fp (9) + alphaP (1)
                           "Unexpected number of state variables for FiniteStrainJ2Plasticity in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  Eigen::VectorXd initialStateVars = Eigen::VectorXd::Zero( nStateVars );
  // Fp = identity (indices 0,4,8 of the flattened 3x3 plastic deformation gradient)
  initialStateVars( 0 ) = initialStateVars( 4 ) = initialStateVars( 8 ) = 1.0;
  initialStateVars( 9 )                                                 = 12.5; // alphaP

  solver.setInitialState( Tensor33d( 0.0 ), initialStateVars );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-6; // tiny: stays elastic, so alphaP is not itself modified
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );

  solver.addStep( step );
  solver.solve();

  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stateVars( 9 ), 12.5, 1e-10 ),
                           "setInitialState() did not seed the solver's initial state variables in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveStep(): retries with a halved time step on StressUpdateFailed, then gives up with
// SolverTimestepExhausted once the time step can no longer be halved below dTMin.
// ---------------------------------------------------------------------------------------------
void testSolveStepRetriesThenThrowsSolverTimestepExhausted()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTALWAYSFAILINGFINITESTRAIN";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-3;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.dTStart                      = 0.1;
  step.dTMin                        = 0.05;

  solver.addStep( step );

  bool threw = false;
  try {
    solver.solve();
  }
  catch ( const Marmot::SolverTimestepExhausted& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "solveStep() must throw SolverTimestepExhausted once the time step cannot be halved "
                           "below dTMin in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveStep(): throws SolverIncrementsExhausted when maxIncrements is reached before timeEnd.
// ---------------------------------------------------------------------------------------------
void testSolveStepThrowsSolverIncrementsExhausted()
{
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-3;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 1.0;
  step.dTStart                      = 0.1; // 10 increments needed to reach timeEnd
  step.maxIncrements                = 2;   // but only 3 are allowed (counter <= maxIncrements)

  solver.addStep( step );

  bool threw = false;
  try {
    solver.solve();
  }
  catch ( const Marmot::SolverIncrementsExhausted& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "solveStep() must throw SolverIncrementsExhausted when maxIncrements is reached "
                           "before timeEnd in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveStep(): the very first attempted increment always uses dT == 0 (a trivial pass, before
// dTStart is assigned on the first successful increment); if dTStart itself would overshoot
// timeEnd, that overshoot must be capped independently of the general (dT-already-set) overshoot
// check earlier in the loop, which ran before dT was reassigned to dTStart.
// ---------------------------------------------------------------------------------------------
void testSolveStepCapsNewlyAssignedDTStartToAvoidOvershoot()
{
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-4;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.timeStart                    = 0.0;
  step.timeEnd                      = 0.05; // shorter than dTStart
  step.dTStart                      = 0.1;  // would overshoot timeEnd if not capped

  solver.addStep( step );
  solver.solve(); // must not throw

  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().time, step.timeEnd, 1e-12 ),
                           "solveStep() did not reach timeEnd exactly after capping the newly-assigned dTStart in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveStep(): a StressUpdateFailed on a genuinely nonzero time increment must halve that
// (already-nonzero) dT and retry, rather than immediately hitting dTMin -- unlike a material that
// fails even on the trivial dT==0 first pass.
// ---------------------------------------------------------------------------------------------
void testSolveStepHalvesDTOnRetryBeforeExhaustingIt()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTFAILSFORNONZERODTFINITESTRAIN";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = Tensor33d( 0.0 );
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.dTStart                     = 0.1;
  step.dTMin                       = 0.05; // exactly one halving (0.1 -> 0.05) before exhaustion

  solver.addStep( step );

  bool threw = false;
  try {
    solver.solve();
  }
  catch ( const Marmot::SolverTimestepExhausted& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "solveStep() must throw SolverTimestepExhausted after halving dT once it is exhausted "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveIncrement(): throws SolverConvergenceFailed with a "NaN encountered" message when the
// tangent is singular.
// ---------------------------------------------------------------------------------------------
void testSolveIncrementThrowsOnNaN()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTSINGULARTANGENTFINITESTRAIN";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget          = Tensor33d( 0.0 );
  step.stressIncrementTarget         = Tensor33d( 0.0 );
  step.stressIncrementTarget( 0, 0 ) = 100.0; // the material always reports tau(0,0) == 5.0
  step.isGradUComponentControlled    = Tensor33t< bool >( false );
  step.isStressComponentControlled   = Tensor33t< bool >( true );

  solver.addStep( step );

  bool threw = false;
  try {
    solver.solve();
  }
  catch ( const Marmot::SolverConvergenceFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "solveIncrement() must throw SolverConvergenceFailed for a singular tangent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveIncrement(): throws SolverConvergenceFailed after exhausting maxIterations when the
// residual can never be driven to zero, even with a well-conditioned (non-singular) tangent.
// ---------------------------------------------------------------------------------------------
void testSolveIncrementThrowsWhenIterationsExhausted()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTNEVERCONVERGINGFINITESTRAINSOLVER";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  solveropts.maxIterations                 = 3; // keep the test fast
  auto solver                              = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget          = Tensor33d( 0.0 );
  step.stressIncrementTarget         = Tensor33d( 0.0 );
  step.stressIncrementTarget( 0, 0 ) = 100.0; // the material always reports tau(0,0) == 5.0
  step.isGradUComponentControlled    = Tensor33t< bool >( false );
  step.isStressComponentControlled   = Tensor33t< bool >( true );

  solver.addStep( step );

  bool threw = false;
  try {
    solver.solve();
  }
  catch ( const Marmot::SolverConvergenceFailed& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "solveIncrement() must throw SolverConvergenceFailed when the residual can never be "
                           "driven to zero in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// printHistory() / exportHistoryToCSV()
// ---------------------------------------------------------------------------------------------
void testPrintHistoryDoesNotThrow()
{
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-3;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );

  solver.addStep( step );
  solver.solve();

  solver.printHistory(); // purely a coverage/smoke check: must not throw
}

void testExportHistoryToCSVWritesExpectedContentWithNoStateVars()
{
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE"; // no state variables
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-3;
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );
  step.dTStart                      = 0.5; // 2 increments

  solver.addStep( step );
  solver.solve();

  const std::string filename = "test_export_history_finitestrain.csv";
  solver.exportHistoryToCSV( filename );

  std::ifstream file( filename );
  throwExceptionOnFailure( file.is_open(),
                           "exportHistoryToCSV() did not create the expected file in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string headerLine;
  std::getline( file, headerLine );
  throwExceptionOnFailure( headerLine.find( "Time" ) != std::string::npos &&
                             headerLine.find( "tau11" ) != std::string::npos &&
                             headerLine.find( "F11" ) != std::string::npos &&
                             headerLine.find( "SV1" ) == std::string::npos,
                           "exportHistoryToCSV() header does not match the expected columns in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  int         dataLines = 0;
  std::string line;
  while ( std::getline( file, line ) )
    if ( !line.empty() )
      dataLines++;

  throwExceptionOnFailure( dataLines == static_cast< int >( solver.getHistory().size() ),
                           "exportHistoryToCSV() did not write one line per history entry in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  file.close();
  std::remove( filename.c_str() );
}

void testExportHistoryToCSVWritesStateVarColumnsForAMaterialWithStateVars()
{
  std::vector< double > materialProperties = { 175000., 80800., 260., 580., 9., 70., 3. }; // purely elastic step
  std::string           matName            = "FINITESTRAINJ2PLASTICITY";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  throwExceptionOnFailure( solver.getNumberOfStateVariables() > 0,
                           "FiniteStrainJ2Plasticity must require at least one state variable in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget         = Tensor33d( 0.0 );
  step.gradUIncrementTarget( 0, 0 ) = 1e-6; // tiny: stays elastic
  step.stressIncrementTarget        = Tensor33d( 0.0 );
  step.isGradUComponentControlled   = Tensor33t< bool >( true );
  step.isStressComponentControlled  = Tensor33t< bool >( false );

  solver.addStep( step );
  solver.solve();

  const std::string filename = "test_export_history_finitestrain_j2.csv";
  solver.exportHistoryToCSV( filename );

  std::ifstream file( filename );
  throwExceptionOnFailure( file.is_open(),
                           "exportHistoryToCSV() did not create the expected file in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string headerLine;
  std::getline( file, headerLine );
  throwExceptionOnFailure( headerLine.find( "SV1" ) != std::string::npos,
                           "exportHistoryToCSV() header is missing the state-variable columns in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string dataLine;
  std::getline( file, dataLine );
  throwExceptionOnFailure( !dataLine.empty(),
                           "exportHistoryToCSV() did not write a data row in " + std::string( __PRETTY_FUNCTION__ ) );

  file.close();
  std::remove( filename.c_str() );
}

void testExportHistoryToCSVThrowsForInvalidPath()
{
  std::vector< double > materialProperties = { 3500., 1500. };
  std::string           matName            = "COMPRESSIBLENEOHOOKE";
  auto                  solveropts         = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverFiniteStrain( matName,
                                                       materialProperties.data(),
                                                       materialProperties.size(),
                                                       solveropts );

  bool threw = false;
  try {
    solver.exportHistoryToCSV( "/nonexistent_directory_xyz/out.csv" );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "exportHistoryToCSV() must throw for a path that cannot be opened for writing in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testResetToInitialStateUndoesAccumulatedDeformation,
    testSetInitialStateSeedsStateVariables,
    testSolveStepRetriesThenThrowsSolverTimestepExhausted,
    testSolveStepThrowsSolverIncrementsExhausted,
    testSolveStepCapsNewlyAssignedDTStartToAvoidOvershoot,
    testSolveStepHalvesDTOnRetryBeforeExhaustingIt,
    testSolveIncrementThrowsOnNaN,
    testSolveIncrementThrowsWhenIterationsExhausted,
    testPrintHistoryDoesNotThrow,
    testExportHistoryToCSVWritesExpectedContentWithNoStateVars,
    testExportHistoryToCSVWritesStateVarColumnsForAMaterialWithStateVars,
    testExportHistoryToCSVThrowsForInvalidPath,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
