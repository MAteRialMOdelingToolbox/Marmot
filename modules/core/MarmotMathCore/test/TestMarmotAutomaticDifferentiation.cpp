#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

// test automatic differentiation for scalar functions
void testAutomaticDifferentiationForScalars()
{
  using namespace Marmot::AutomaticDifferentiation;

  // check scalar first derivative
  auto scalar_func = [&]( const autodiff::dual x ) {
    const autodiff::dual res = x * x * exp( x );
    return res;
  };

  const double x_     = 2.5;
  const double df_dx_ = 2. * x_ * exp( x_ ) + x_ * x_ * exp( x_ );
  // check scalar first derivative using 1st order duals
  throwExceptionOnFailure( checkIfEqual( df_dx( scalar_func, x_ ), df_dx_ ), "Error in df_dx with autodiff::dual" );

  // check scalar first derivative
  auto scalar_func_2nd = [&]( const autodiff::dual2nd x ) {
    const autodiff::dual2nd res = x * x * exp( x );
    return res;
  };

  // check scalar first derivative using 2nd order duals
  const autodiff::dual xDual     = 2.5;
  const autodiff::dual df_dx_2nd = 2. * x_ * exp( xDual ) + xDual * xDual * exp( xDual );
  throwExceptionOnFailure( checkIfEqual( df_dx( scalar_func_2nd, xDual ), df_dx_2nd ),
                           "Error in df_dx with autodiff::dual2nd" );
}

void testADForVectorValuedFunctions()
{
  using namespace Marmot::AutomaticDifferentiation;

  // check jacobian with autodiff::dual
  auto vector_func = [&]( const autodiff::VectorXdual x ) {
    const autodiff::VectorXdual res = 2 * x;
    return res;
  };

  VectorXd X( 3 );
  X << 1, 2, 3;

  auto [F, J] = dF_dX( vector_func, X );

  throwExceptionOnFailure( checkIfEqual< double >( F, 2 * X ),
                           "Error in vector function evaluation for jacobian computation with autodiff::dual" );
  throwExceptionOnFailure( checkIfEqual< double >( J, 2 * MatrixXd::Identity( 3, 3 ) ),
                           "Error in vector function jacobian with autodiff::dual" );

  // check jacobian with autodiff::dual2nd
  auto vector_func_2nd = [&]( const autodiff::VectorXdual2nd x ) {
    const autodiff::VectorXdual2nd res = 2 * x;
    return res;
  };

  VectorXdual XDual( 3 );
  XDual( 0 ).val  = 1;
  XDual( 0 ).grad = 1;
  XDual( 1 ).val  = 2;
  XDual( 1 ).grad = 1;
  XDual( 2 ).val  = 3;
  XDual( 2 ).grad = 1;

  auto [F_, J_] = dF_dX_2nd( vector_func_2nd, XDual );

  MatrixXdual target = 2 * MatrixXdual::Identity( 3, 3 );

  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( F_, 2. * XDual ),
                           "Error in vector function evaluation for jacobian computation with autodiff::dual2nd" );
  throwExceptionOnFailure( checkIfEqual< autodiff::dual >( J, target ),
                           "Error in vector function jacobian with autodiff::dual2nd" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testAutomaticDifferentiationForScalars,
                                                       testADForVectorValuedFunctions };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
