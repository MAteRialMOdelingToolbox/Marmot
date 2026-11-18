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

  // check jacobian with dual2nd: F( x ) = ( x0 x1, x1^2, x0 + x2^3 ), with x carrying its derivative with respect to an
  // outer variable t (dx/dt = s); F_ and J_ must carry dF/dt and dJ/dt
  auto vector_func_2nd = [&]( const VectorXdual2nd& x ) {
    VectorXdual2nd res( 3 );
    res( 0 ) = x( 0 ) * x( 1 );
    res( 1 ) = x( 1 ) * x( 1 );
    res( 2 ) = x( 0 ) + x( 2 ) * x( 2 ) * x( 2 );
    return res;
  };

  const double x[3] = { 1.5, -0.5, 2.0 }, s[3] = { 0.5, -1.0, 2.0 };
  VectorXdual  XDual( 3 );
  for ( int i = 0; i < 3; i++ ) {
    XDual( i ) = dual( x[i] );
    seed< 1 >( XDual( i ), s[i] );
  }

  auto [F_, J_] = dF_dX_2nd( vector_func_2nd, XDual );

  const double F2[3]     = { x[0] * x[1], x[1] * x[1], x[0] + x[2] * x[2] * x[2] };
  const double dF2[3]    = { s[0] * x[1] + x[0] * s[1], 2 * x[1] * s[1], s[0] + 3 * x[2] * x[2] * s[2] };
  const double J2[3][3]  = { { x[1], x[0], 0 }, { 0, 2 * x[1], 0 }, { 1, 0, 3 * x[2] * x[2] } };
  const double dJ2[3][3] = { { s[1], s[0], 0 }, { 0, 2 * s[1], 0 }, { 0, 0, 6 * x[2] * s[2] } };
  for ( int i = 0; i < 3; i++ ) {
    throwExceptionOnFailure( checkIfEqual( value( F_( i ) ), F2[i], 1e-14 ) &&
                               checkIfEqual( derivative< 1 >( F_( i ) ), dF2[i], 1e-14 ),
                             "Error in vector function evaluation for jacobian computation with dual2nd" );
    for ( int j = 0; j < 3; j++ )
      throwExceptionOnFailure( checkIfEqual( value( J_( i, j ) ), J2[i][j], 1e-14 ) &&
                                 checkIfEqual( derivative< 1 >( J_( i, j ) ), dJ2[i][j], 1e-14 ),
                               "Error in vector function jacobian with dual2nd" );
  }
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
