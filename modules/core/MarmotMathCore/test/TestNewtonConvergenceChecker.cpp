#include "Marmot/MarmotTesting.h"
#include "Marmot/NewtonConvergenceChecker.h"
#include <Eigen/Core>

using namespace Marmot::Testing;

int main()
{

  int    nMaxNewtonCycles    = 10;
  int    nMaxNewtonCyclesAlt = 20;
  double newtonTol           = 1e-6;
  double newtonRTol          = 1e-6;
  double newtonTolAlt        = 1e-3;
  double newtonRTolAlt       = 1e-3;

  Eigen::VectorXd residualScaleVector( 1 );
  residualScaleVector << 1.;

  Marmot::NumericalAlgorithms::NewtonConvergenceChecker checker( residualScaleVector,
                                                                 nMaxNewtonCycles,
                                                                 nMaxNewtonCyclesAlt,
                                                                 newtonTol,
                                                                 newtonRTol,
                                                                 newtonTolAlt,
                                                                 newtonRTolAlt );

  // the "last newton increment"
  Eigen::VectorXd increment( 1 );
  increment << 1.;

  // the "solution vector"
  Eigen::VectorXd X( 1 );
  X << 1.;

  throwExceptionOnFailure( checkIfEqual( checker.relativeNorm( increment, X ), 1.0 ),
                           "checker.relativeNorm( increment, reference ) != 1." );

  // the "residual vector"
  Eigen::VectorXd residual( 1 );
  residual << 1.;

  throwExceptionOnFailure( checkIfEqual( checker.residualNorm( residual ), 1. ),
                           "checker.residualNorm( residual ) != 1." );

  residual << 1.;
  throwExceptionOnFailure( checker.iterationFinished( residual, X, increment, 5 ) == false,
                           "checker.iterationFinished( residual, X, increment, 0 ) != false" );

  throwExceptionOnFailure( checker.isConverged( residual, X, increment, 5 ) == false,
                           "checker.isConverged( residual, X, increment, 0 ) != false" );

  // the "residual vector"
  residual << 1e-5;
  increment << 1e-12;

  // with normal tol it should not pass
  throwExceptionOnFailure( checker.iterationFinished( residual, X, increment, 5 ) == false,
                           "checker.iterationFinished( residual, X, increment, 5 ) != false" );

  throwExceptionOnFailure( checker.isConverged( residual, X, increment, 5 ) == false,
                           "checker.isConverged( residual, X, increment, 5 ) != false" );

  // with alt tol it should pass
  throwExceptionOnFailure( checker.iterationFinished( residual, X, increment, 15 ) == true,
                           "checker.iterationFinished( residual, X, increment, 15 ) != true" );

  throwExceptionOnFailure( checker.isConverged( residual, X, increment, 15 ) == true,
                           "checker.isConverged( residual, X, increment, 15 ) != true" );

  // convergence within nMaxNewtonCycles using the *normal* (non-alt) tolerances
  residual << 1e-9;
  increment << 1e-9;
  throwExceptionOnFailure( checker.isConverged( residual, X, increment, 5 ) == true,
                           "checker.isConverged() must be true for a residual/increment satisfying the normal "
                           "tolerances within nMaxNewtonCycles" );
  throwExceptionOnFailure( checker.iterationFinished( residual, X, increment, 5 ) == true,
                           "checker.iterationFinished() must be true once isConverged() is true" );

  // relativeNorm: increment norm below 1e-14 must short-circuit to incNorm itself
  Eigen::VectorXd tinyIncrement( 1 );
  tinyIncrement << 1e-15;
  throwExceptionOnFailure( checkIfEqual( checker.relativeNorm( tinyIncrement, X ), tinyIncrement.norm() ),
                           "relativeNorm() must return incNorm unchanged when incNorm < 1e-14" );

  // relativeNorm: reference norm below 1e-12 must yield 0.0 (undefined relative norm)
  Eigen::VectorXd tinyX( 1 );
  tinyX << 1e-13;
  throwExceptionOnFailure( checkIfEqual( checker.relativeNorm( increment, tinyX ), 0.0 ),
                           "relativeNorm() must return 0.0 when the reference norm < 1e-12" );

  return 0;
}
