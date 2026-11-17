#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

// test automatic differentiation for scalar functions
void testAutomaticDifferentiationForScalars()
{
  using namespace Marmot::Math::Differentiation::Automatic;

  // check scalar first derivative
  auto scalar_func = [&]( const dual x ) {
    const dual res = x * x * exp( x );
    return res;
  };

  const double x_     = 2.5;
  const double df_dx_ = 2. * x_ * exp( x_ ) + x_ * x_ * exp( x_ );
  // check scalar first derivative using 1st order duals
  throwExceptionOnFailure( checkIfEqual( df_dx( scalar_func, x_ ), df_dx_ ), "Error in df_dx with dual" );

  // check scalar first derivative
  auto scalar_func_2nd = [&]( const dual2nd x ) {
    const dual2nd res = x * x * exp( x );
    return res;
  };

  // check scalar first derivative using 2nd order duals
  const dual xDual     = 2.5;
  const dual df_dx_2nd = 2. * x_ * exp( xDual ) + xDual * xDual * exp( xDual );
  throwExceptionOnFailure( checkIfEqual( df_dx( scalar_func_2nd, xDual ), df_dx_2nd ), "Error in df_dx with dual2nd" );
}

void testADForVectorValuedFunctions()
{
  using namespace Marmot::Math::Differentiation::Automatic;

  // check jacobian with dual
  auto vector_func = [&]( const VectorXdual x ) {
    const VectorXdual res = 2 * x;
    return res;
  };

  Eigen::VectorXd X( 3 );
  X << 1, 2, 3;

  auto [F, J] = dF_dX( vector_func, X );

  throwExceptionOnFailure( checkIfEqual< double >( F, 2 * X ),
                           "Error in vector function evaluation for jacobian computation with dual" );
  throwExceptionOnFailure( checkIfEqual< double >( J, 2 * Eigen::MatrixXd::Identity( 3, 3 ) ),
                           "Error in vector function jacobian with dual" );

  // check jacobian with dual2nd
  auto vector_func_2nd = [&]( const VectorXdual2nd x ) {
    const VectorXdual2nd res = 2 * x;
    return res;
  };

  VectorXdual XDual( 3 );
  XDual( 0 ) = dual( 1 );
  seed< 1 >( XDual( 0 ), 1 );
  XDual( 1 ) = dual( 2 );
  seed< 1 >( XDual( 1 ), 1 );
  XDual( 2 ) = dual( 3 );
  seed< 1 >( XDual( 2 ), 1 );

  auto [F_, J_] = dF_dX_2nd( vector_func_2nd, XDual );

  MatrixXdual target = 2 * MatrixXdual::Identity( 3, 3 );

  throwExceptionOnFailure( checkIfEqual< dual >( F_, 2. * XDual ),
                           "Error in vector function evaluation for jacobian computation with dual2nd" );
  throwExceptionOnFailure( checkIfEqual< dual >( J, target ), "Error in vector function jacobian with dual2nd" );
}

void testShiftTo2ndOrderDualForVectors()
{
  using namespace Marmot::Math::Differentiation::Automatic;

  // seed both the value AND the first-order derivative of every input, so that
  // this test actually exercises shiftTo2ndOrderDual's documented behavior of
  // moving each input's first-order derivative into dual2nd.val.grad. With an
  // all-zero seed (as before), a broken implementation that drops the incoming
  // gradient entirely would still pass this test.
  VectorXdual X( 3 );
  X( 0 ) = dual( 1.5 );
  seed< 1 >( X( 0 ), 0.5 );
  X( 1 ) = dual( -2.5 );
  seed< 1 >( X( 1 ), -1.25 );
  X( 2 ) = dual( 3.5 );
  seed< 1 >( X( 2 ), 2.0 );

  const VectorXdual2nd X2nd = shiftTo2ndOrderDual( X );

  throwExceptionOnFailure( static_cast< int >( X2nd.size() ) == 3,
                           "shiftTo2ndOrderDual(VectorXdual) returned the wrong size in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  for ( int i = 0; i < 3; i++ ) {
    // ground truth: the already-tested scalar overload of shiftTo2ndOrderDual
    const dual2nd expected = shiftTo2ndOrderDual( X( i ) );

    throwExceptionOnFailure( checkIfEqual( value( X2nd( i ) ), value( expected ) ),
                             "shiftTo2ndOrderDual(VectorXdual) does not match the scalar overload's value at index " +
                               std::to_string( i ) + " in " + std::string( __PRETTY_FUNCTION__ ) );

    throwExceptionOnFailure( checkIfEqual( coefficients< 2 >( X2nd( i ) )[1], coefficients< 2 >( expected )[1] ),
                             "shiftTo2ndOrderDual(VectorXdual) does not match the scalar overload's "
                             "first-order derivative (val.grad) at index " +
                               std::to_string( i ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testAutomaticDifferentiationForScalars,
                                                       testADForVectorValuedFunctions,
                                                       testShiftTo2ndOrderDualForVectors };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
