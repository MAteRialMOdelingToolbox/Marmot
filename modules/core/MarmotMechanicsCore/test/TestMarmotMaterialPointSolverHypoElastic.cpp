#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Solvers;

namespace {

  // A material that always fails to update its stress, used to exercise solveStep()'s
  // cutback-retry-then-give-up logic deterministically (no real material fails unconditionally).
  class AlwaysFailingHypoElasticMaterial : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    void computeStress( state3D&, Marmot::Matrix6d&, const Marmot::Vector6d&, const timeInfo& ) const override
    {
      throw Marmot::StressUpdateFailed( "test-only material always fails" );
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  // A material whose stress output is a strain-independent constant with a zero tangent, used to
  // exercise solveIncrement()'s own Newton-loop non-convergence throw: for a stress-controlled
  // component, the residual can then never be driven to zero, regardless of how many correction
  // steps are attempted.
  class NeverConvergingSolverMaterial : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    void computeStress( state3D&          state,
                        Marmot::Matrix6d& dStress_dStrain,
                        const Marmot::Vector6d&,
                        const timeInfo& ) const override
    {
      state.stress      = Marmot::Vector6d::Zero();
      state.stress( 0 ) = 5.0; // never equal to any nonzero target, and independent of dStrain
      dStress_dStrain.setZero();
    }

    double getDensity( const double* ) const override { return 1.0; }
  };

  const bool alwaysFailingRegistered = MarmotLibrary::MarmotMaterialHypoElasticFactory::registerMaterial<
    AlwaysFailingHypoElasticMaterial >( "TESTALWAYSFAILINGHYPOELASTIC" );
  const bool neverConvergingRegistered = MarmotLibrary::MarmotMaterialHypoElasticFactory::registerMaterial<
    NeverConvergingSolverMaterial >( "TESTNEVERCONVERGINGHYPOELASTICSOLVER" );

} // namespace

// ---------------------------------------------------------------------------------------------
// setInitialState() / resetToInitialState()
// ---------------------------------------------------------------------------------------------
void testSetInitialStateAndResetToInitialState()
{
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  Marmot::Vector6d initialStress = Marmot::Vector6d::Zero();
  initialStress( 0 )             = 42.0;
  Eigen::VectorXd initialStateVars( solver.getNumberOfStateVariables() );
  initialStateVars.setZero();

  solver.setInitialState( initialStress, initialStateVars );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget       = Marmot::Vector6d::Zero(); // zero strain increment: stress must stay at 42

  solver.addStep( step );
  solver.solve();

  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress( 0 ), 42.0, 1e-10 ),
                           "setInitialState() did not seed the solver's starting stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  solver.resetToInitialState();
  throwExceptionOnFailure( solver.getHistory().empty(),
                           "resetToInitialState() must clear the recorded history in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // After resetting, solving the same step again must reproduce the same starting stress.
  // resetToInitialState() does not clear queued steps, so avoid re-solving the (already-queued)
  // first step on top of itself.
  solver.clearSteps();
  solver.addStep( step );
  solver.solve();
  throwExceptionOnFailure( checkIfEqual( solver.getHistory().back().stress( 0 ), 42.0, 1e-10 ),
                           "resetToInitialState() did not restore the initial stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveStep(): retries with a halved time step on StressUpdateFailed, then gives up with
// SolverTimestepExhausted once the time step can no longer be halved below dTMin.
// ---------------------------------------------------------------------------------------------
void testSolveStepRetriesThenThrowsSolverTimestepExhausted()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTALWAYSFAILINGHYPOELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.dTStart                     = 0.1;
  step.dTMin                       = 0.05;

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
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.timeStart                   = 0.0;
  step.timeEnd                     = 1.0;
  step.dTStart                     = 0.1; // 10 increments needed to reach timeEnd
  step.maxIncrements               = 2;   // but only 3 are allowed (counter <= maxIncrements)

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
// solveIncrement(): a mixed strain/stress-controlled increment that is not already exactly
// satisfied by the initial guess needs at least one Newton correction step -- unlike every
// fully-strain-controlled LinearElastic test elsewhere, which converges in "0 iterations" because
// the initial guess already is the exact answer.
// ---------------------------------------------------------------------------------------------
void testSolveIncrementMixedControlNeedsNewtonCorrection()
{
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { false, true, true, true, true, true };
  step.isStressComponentControlled = { true, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.stressIncrementTarget( 0 )  = 240.0;

  solver.addStep( step );
  solver.solve();

  const auto  history = solver.getHistory();
  const auto& last    = history.back();

  throwExceptionOnFailure( checkIfEqual( last.stress( 0 ), 240.0, 1e-6 ),
                           "solveIncrement() did not converge to the target stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Cross-check against a freshly constructed material evaluated directly at the converged
  // strain, rather than re-deriving the expected stress analytically.
  Marmot::Materials::LinearElastic    mat( materialProperties.data(), materialProperties.size(), 1 );
  MarmotMaterialHypoElastic&          matBase = mat;
  MarmotMaterialHypoElastic::state3D  state( Marmot::Vector6d::Zero(), 0.0, 0.0, nullptr );
  Marmot::Matrix6d                    tangent;
  MarmotMaterialHypoElastic::timeInfo timeInfo{ 0.0, 1.0 };
  matBase.computeStress( state, tangent, last.strain, timeInfo );

  throwExceptionOnFailure( checkIfEqual< double >( state.stress, last.stress, 1e-6 ),
                           "The solver's converged (stress, strain) pair is not self-consistent with a direct "
                           "material evaluation in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// solveIncrement(): throws SolverConvergenceFailed when the Newton loop cannot drive the residual
// below tolerance within maxIterations.
// ---------------------------------------------------------------------------------------------
void testSolveIncrementThrowsSolverConvergenceFailed()
{
  std::vector< double > materialProperties = { 0.0 };
  std::string           matName            = "TESTNEVERCONVERGINGHYPOELASTICSOLVER";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  solveropts.maxIterations                 = 3; // keep the test fast
  auto solver                              = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { false, true, true, true, true, true };
  step.isStressComponentControlled = { true, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.stressIncrementTarget( 0 )  = 100.0; // the material always reports stress(0) == 5.0

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
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();

  solver.addStep( step );
  solver.solve();

  solver.printHistory(); // purely a coverage/smoke check: must not throw
}

void testExportHistoryToCSVWritesExpectedContent()
{
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();
  step.dTStart                     = 0.5; // 2 increments

  solver.addStep( step );
  solver.solve();

  const std::string filename = "test_export_history_hypoelastic.csv";
  solver.exportHistoryToCSV( filename );

  std::ifstream file( filename );
  throwExceptionOnFailure( file.is_open(),
                           "exportHistoryToCSV() did not create the expected file in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string headerLine;
  std::getline( file, headerLine );
  throwExceptionOnFailure( headerLine.find( "Time" ) != std::string::npos &&
                             headerLine.find( "Stress_11" ) != std::string::npos &&
                             headerLine.find( "StateVar_1" ) == std::string::npos, // LinearElastic has none
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

void testExportHistoryToCSVThrowsForInvalidPath()
{
  std::vector< double > materialProperties = { 20000., 0.25 };
  std::string           matName            = "LINEARELASTIC";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
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

// ---------------------------------------------------------------------------------------------
// exportHistoryToCSV(): the state-variable columns (both header and per-row) are only written for
// a material that actually has state variables -- LinearElastic (used above) has none.
// ---------------------------------------------------------------------------------------------
void testExportHistoryToCSVWritesStateVarColumnsForAMaterialWithStateVars()
{
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20. }; // VonMises, purely elastic
  std::string           matName            = "VONMISES";
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                  solver             = MarmotMaterialPointSolverHypoElastic( matName,
                                                      materialProperties.data(),
                                                      materialProperties.size(),
                                                      solveropts );

  MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-4; // small: stays elastic, well below the yield stress
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();

  solver.addStep( step );
  solver.solve();

  throwExceptionOnFailure( solver.getNumberOfStateVariables() > 0,
                           "VonMises must require at least one state variable in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const std::string filename = "test_export_history_vonmises.csv";
  solver.exportHistoryToCSV( filename );

  std::ifstream file( filename );
  throwExceptionOnFailure( file.is_open(),
                           "exportHistoryToCSV() did not create the expected file in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string headerLine;
  std::getline( file, headerLine );
  throwExceptionOnFailure( headerLine.find( "StateVar_1" ) != std::string::npos,
                           "exportHistoryToCSV() header is missing the state-variable columns in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::string dataLine;
  std::getline( file, dataLine );
  throwExceptionOnFailure( !dataLine.empty(),
                           "exportHistoryToCSV() did not write a data row in " + std::string( __PRETTY_FUNCTION__ ) );

  file.close();
  std::remove( filename.c_str() );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testSetInitialStateAndResetToInitialState,
    testSolveStepRetriesThenThrowsSolverTimestepExhausted,
    testSolveStepThrowsSolverIncrementsExhausted,
    testSolveIncrementMixedControlNeedsNewtonCorrection,
    testSolveIncrementThrowsSolverConvergenceFailed,
    testPrintHistoryDoesNotThrow,
    testExportHistoryToCSVWritesExpectedContent,
    testExportHistoryToCSVWritesStateVarColumnsForAMaterialWithStateVars,
    testExportHistoryToCSVThrowsForInvalidPath,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
