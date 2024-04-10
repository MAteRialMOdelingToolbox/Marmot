#include "Marmot/MarmotNumericalIntegration.h"
#include <iostream>
#include <iterator>

using namespace Marmot::NumericalAlgorithms::Integration;

int main()
{

  auto f = [&]( double x ) { return x * x * x; };
  
  // Test Simpson Rule
  double res = integrateScalarFunction( f, { 0., 1. }, 1, integrationRule::simpson );

  if ( res != 0.25 ) {
    std::cout << "Numerical Integration with Simpson Rule failed: " << res << " != 0.25" << std::endl;
    return 1;
  }

  return 0;
}
