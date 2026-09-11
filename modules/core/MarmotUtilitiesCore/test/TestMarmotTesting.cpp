#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <limits>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

// ---------------------------------------------------------------------------------------------
// Most of Marmot::Testing's own source lines are only reached when an equality check actually
// *fails* -- which, correctly, never happens in a clean, all-passing suite. Testing them means
// deliberately triggering a "failure" and checking the framework reports/returns it correctly,
// without letting that propagate as a real failure of this test suite.
// ---------------------------------------------------------------------------------------------

void testGetStringFormatsDoubleAndDual()
{
  throwExceptionOnFailure( getString( 1.5 ) == std::to_string( 1.5 ),
                           "getString(double) failed in " + std::string( __PRETTY_FUNCTION__ ) );

  autodiff::dual d    = 2.0;
  d.grad              = 3.0;
  const std::string s = getString( d );
  throwExceptionOnFailure( s.find( std::to_string( 2.0 ) ) != std::string::npos &&
                             s.find( std::to_string( 3.0 ) ) != std::string::npos,
                           "getString(autodiff::dual) failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

void testCheckIfEqualReturnsFalseForNaNAndInf()
{
  const double nan = std::numeric_limits< double >::quiet_NaN();
  const double inf = std::numeric_limits< double >::infinity();

  throwExceptionOnFailure( checkIfEqual( nan, 1.0 ) == false,
                           "checkIfEqual() must return false when the first argument is NaN in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( 1.0, nan ) == false,
                           "checkIfEqual() must return false when the second argument is NaN in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( inf, 1.0 ) == false,
                           "checkIfEqual() must return false when the first argument is inf in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( 1.0, inf ) == false,
                           "checkIfEqual() must return false when the second argument is inf in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testCheckIfEqualReturnsFalseForMismatchedScalars()
{
  throwExceptionOnFailure( checkIfEqual( 1.0, 2.0, 1e-10 ) == false,
                           "checkIfEqual(double) must return false for clearly different values in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const autodiff::dual a = 1.0, b = 2.0;
  throwExceptionOnFailure( checkIfEqual( a, b, 1e-10 ) == false,
                           "checkIfEqual(autodiff::dual) must return false for different values in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( a, a, 1e-10 ) == true,
                           "checkIfEqual(autodiff::dual) must return true for identical values in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  const std::complex< double > c1( 1.0, 0.0 ), c2( 1.0, 1.0 );
  throwExceptionOnFailure( checkIfEqual( c1, c2, 1e-10 ) == false,
                           "checkIfEqual(complex<double>) must return false for different imaginary parts in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( c1, c1, 1e-10 ) == true,
                           "checkIfEqual(complex<double>) must return true for identical values in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testCheckIfEqualMatrixHintPathForMismatch()
{
  // Exercises the "-> HINT: a(i,j) = ... != b(i,j) = ..." diagnostic path (and, through it,
  // getString()), which is only reached inside a *failing* matrix comparison.
  Eigen::MatrixXd a( 1, 1 ), b( 1, 1 );
  a( 0, 0 ) = 1.0;
  b( 0, 0 ) = 2.0;

  throwExceptionOnFailure( checkIfEqual< double >( a, b, 1e-10 ) == false,
                           "checkIfEqual() for mismatched matrices must return false in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testThrowExceptionOnFailureThrowsForFalseCondition()
{
  bool threw = false;
  try {
    throwExceptionOnFailure( false, "deliberate test failure message" );
  }
  catch ( const std::runtime_error& e ) {
    threw = std::string( e.what() ) == "deliberate test failure message";
  }
  throwExceptionOnFailure( threw,
                           "throwExceptionOnFailure(false, ...) must throw std::runtime_error with the given "
                           "message in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Must not throw for a true condition.
  throwExceptionOnFailure( true, "must not be thrown" );
}

void testExecuteTestsAndCollectExceptionsAggregatesFailures()
{
  bool passed = false;

  const std::vector< std::function< void() > > innerTests = {
    []() { /* passes */ },
    []() { throw std::runtime_error( "deliberate inner failure" ); },
  };

  bool threw = false;
  try {
    executeTestsAndCollectExceptions( innerTests );
  }
  catch ( const std::runtime_error& e ) {
    threw = std::string( e.what() ) == "some tests failed";
  }
  throwExceptionOnFailure( threw,
                           "executeTestsAndCollectExceptions() must throw \"some tests failed\" when a test "
                           "throws in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // A fully passing set must not throw.
  const std::vector< std::function< void() > > passingTests = { []() {}, [&passed]() { passed = true; } };
  executeTestsAndCollectExceptions( passingTests );
  throwExceptionOnFailure( passed,
                           "executeTestsAndCollectExceptions() did not run a passing test in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

namespace {

  Marmot::Solvers::MarmotMaterialPointSolverHypoElastic makeIsotropicSolver()
  {
    static std::vector< double > materialProperties = { 20000., 0.25 };
    static std::string           matName            = "LINEARELASTIC";
    auto                         solveropts = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::SolverOptions();
    return Marmot::Solvers::MarmotMaterialPointSolverHypoElastic( matName,
                                                                  materialProperties.data(),
                                                                  materialProperties.size(),
                                                                  solveropts );
  }

} // namespace

void testSpinTurbokreiselReturnsFalseForNoSteps()
{
  auto solver = makeIsotropicSolver();
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-8, 1e-8 ) == false,
                           "spinTurbokreisel() must return false when no steps were added in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSpinTurbokreiselReturnsFalseForAStressControlledStep()
{
  auto solver = makeIsotropicSolver();

  Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, false };
  step.isStressComponentControlled = { false, false, false, false, false, true }; // not pure strain control
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();

  solver.addStep( step );

  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-8, 1e-8 ) == false,
                           "spinTurbokreisel() must return false for a step that is not purely strain-controlled "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSpinTurbokreiselDetectsAnAnisotropicMaterial()
{
  // A transversely isotropic (i.e. genuinely NOT isotropic) LinearElastic material: E1, E2, nu12,
  // nu23, G12, and the two direction vectors defining the local coordinate system.
  static std::vector< double > materialProperties = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1 };
  static std::string           matName            = "LINEARELASTIC";
  auto                         solveropts = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                         solver     = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic( matName,
                                                                       materialProperties.data(),
                                                                       materialProperties.size(),
                                                                       solveropts );

  Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();

  solver.addStep( step );

  // A transversely isotropic material's response depends on orientation relative to its material
  // axes, so rotating the loading direction (as spinTurbokreisel does) must be detected as a
  // stress and/or tangent mismatch against the un-rotated reference response.
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-6, 1e-6 ) == false,
                           "spinTurbokreisel() must detect a transversely isotropic (non-isotropic) material as "
                           "failing the isotropy check in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSpinTurbokreiselDetectsATangentOnlyMismatch()
{
  // Same anisotropic setup as above, but with a deliberately loose stress tolerance: exercises the
  // separate tangent-mismatch failure branch (reached only once the stress check has already
  // passed), rather than the stress-mismatch branch above.
  static std::vector< double > materialProperties = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1 };
  static std::string           matName            = "LINEARELASTIC";
  auto                         solveropts = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::SolverOptions();
  auto                         solver     = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic( matName,
                                                                       materialProperties.data(),
                                                                       materialProperties.size(),
                                                                       solveropts );

  Marmot::Solvers::MarmotMaterialPointSolverHypoElastic::Step step;
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.strainIncrementTarget       = Marmot::Vector6d::Zero();
  step.strainIncrementTarget( 0 )  = 1e-3;
  step.stressIncrementTarget       = Marmot::Vector6d::Zero();

  solver.addStep( step );

  throwExceptionOnFailure( spinTurbokreisel( solver, 1e2, 1e-6 ) == false,
                           "spinTurbokreisel() must detect a tangent-only mismatch under a loose stress "
                           "tolerance in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testGetStringFormatsDoubleAndDual,
    testCheckIfEqualReturnsFalseForNaNAndInf,
    testCheckIfEqualReturnsFalseForMismatchedScalars,
    testCheckIfEqualMatrixHintPathForMismatch,
    testThrowExceptionOnFailureThrowsForFalseCondition,
    testExecuteTestsAndCollectExceptionsAggregatesFailures,
    testSpinTurbokreiselReturnsFalseForNoSteps,
    testSpinTurbokreiselReturnsFalseForAStressControlledStep,
    testSpinTurbokreiselDetectsAnAnisotropicMaterial,
    testSpinTurbokreiselDetectsATangentOnlyMismatch,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
