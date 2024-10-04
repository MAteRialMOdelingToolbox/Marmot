#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

// test automatic differentiation
void testAutomaticDifferentiationForScalars()
{
  using namespace Marmot::AutomaticDifferentiation;

  // check scalar first derivative
  auto scalar_func = [&]( const autodiff::dual x ) {
    const autodiff::dual res = x * x * exp( x );
    return res;
  };

  const double x_ = 2.5;
  checkIfEqual( df_dx( scalar_func, x_ ), 2. * x_ * exp( x_ ) + x_ * x_ * exp( x_ ) );

  // check scalar first derivative
  auto scalar_func_2nd = [&]( const autodiff::dual2nd x ) {
    const autodiff::dual2nd res = x * x * exp( x );
    return res;
  };

  // check scalar first derivative using 2nd order duals
  const autodiff::dual xDual = 2.5;
  checkIfEqual( df_dx( scalar_func_2nd, xDual ), 2. * x_ * exp( xDual ) + xDual * xDual * exp( xDual ) );
}

int main()
{
  testAutomaticDifferentiationForScalars();
  return 0;
}
