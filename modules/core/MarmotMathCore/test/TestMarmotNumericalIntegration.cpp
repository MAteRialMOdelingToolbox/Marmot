#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <tuple>

using namespace Marmot::NumericalAlgorithms::Integration;
using namespace Marmot::Testing; // Added for testing utilities

void testIntegrateScalarFunction()
{
  /*
   * Test the integration of a scalar function using different integration rules.
   * The function to be integrated is f(x) = x^2, and the integration limits are [0, 1].
   * The expected result for all integration rules is 1/3.
   */
  // Define a simple function to integrate: f(x) = x^2
  auto testFunction = []( const double x ) { return std::pow( x, 2 ); };

  // Define integration limits
  std::tuple< double, double > limits = { 0.0, 1.0 };

  // Define number of integration steps
  int nSteps = 1000;

  // Analytical result: integral of x^2 from 0 to 1 is 1/3
  double analyticalResult = 1.0 / 3.0;
  double tolerance        = 1e-6; // Tolerance for comparison

  // Test Midpoint Rule
  double midpointResult = integrateScalarFunction( testFunction, limits, nSteps, integrationRule::midpoint );
  throwExceptionOnFailure( checkIfEqual( midpointResult, analyticalResult, tolerance ),
                           "Midpoint rule integration failed." );

  // Test Trapezoidal Rule
  double trapezoidalResult = integrateScalarFunction( testFunction, limits, nSteps, integrationRule::trapezodial );
  throwExceptionOnFailure( checkIfEqual( trapezoidalResult, analyticalResult, tolerance ),
                           "Trapezoidal rule integration failed." );

  // Test Simpson's Rule
  double simpsonResult = integrateScalarFunction( testFunction, limits, nSteps, integrationRule::simpson );
  throwExceptionOnFailure( checkIfEqual( simpsonResult, analyticalResult, tolerance ),
                           "Simpson's rule integration failed." );
}

int main()
{
  testIntegrateScalarFunction();

  return 0;
}
