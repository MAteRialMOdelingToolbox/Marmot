#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"

#include "Marmot/MarmotConstants.h"

using namespace Marmot::Testing;
using namespace Marmot::NumericalAlgorithms::Differentiation;

// test numeric differentiation for scalar functions
void testNumericDifferentiationForScalars()
{
  // check scalar first derivative
  auto scalar_func = []( const double x ) {
    const double res = 3.4 * x * x * exp( x ) + 10. * std::cos( x );
    return res;
  };

  const double x_     = 2.5;
  const double df_dx_ = 3.4 * 2. * x_ * exp( x_ ) + 3.4 * x_ * x_ * exp( x_ ) - 10. * std::sin( x_ );

  throwExceptionOnFailure( checkIfEqual( forwardDifference( scalar_func, x_ ), df_dx_, 5e-5 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: forward difference doesn't yield the right result" );
  throwExceptionOnFailure( checkIfEqual( centralDifference( scalar_func, x_ ), df_dx_, 1e-7 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: central difference doesn't yield the right result" );
}

void testNumericalDifferntiationForVectorValuedFunctions()
{
  using namespace Eigen;

  // check jacobian
  auto vector_func = []( const VectorXd x ) {
    const VectorXd res = 2. * ( x.array().pow( 3.0 ) ).matrix() + 10. * x.array().cos().matrix();
    return res;
  };

  VectorXd X( 3 );
  X << 1, 2, 3;

  auto J_forward = forwardDifference( vector_func, X );
  auto J_central = centralDifference( vector_func, X );
  auto df_dx_    = 2. * 3. * X.array().pow( 2.0 ).matrix() + 10. * ( -X.array().sin().matrix() );

  Matrix3d target = Matrix3d::Zero();
  target( 0, 0 )  = df_dx_( 0 );
  target( 1, 1 )  = df_dx_( 1 );
  target( 2, 2 )  = df_dx_( 2 );

  throwExceptionOnFailure( checkIfEqual< double >( J_forward, target, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with forward difference doesn't yield the right result" );
  throwExceptionOnFailure( checkIfEqual< double >( J_central, target, 1e-9 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with central difference doesn't yield the right result" );
}

void testComplexStepDerivativeForScalars()
{
  // check scalar first derivative
  auto scalar_func = []( const std::complex< double > x ) {
    const std::complex< double > res = 3.4 * x * x * exp( x ) + 10. * std::cos( x );
    return res;
  };

  const double x_     = 2.5;
  const double df_dx_ = 3.4 * 2. * x_ * exp( x_ ) + 3.4 * x_ * x_ * exp( x_ ) - 10. * std::sin( x_ );

  throwExceptionOnFailure( checkIfEqual( Complex::forwardDifference( scalar_func, x_ ), df_dx_, 1e-13 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: complex forward difference doesn't yield the right result" );
}

void testComplexStepDerivativeForVectorValuedFunctions()
{
  using namespace Eigen;

  // check jacobian
  auto vector_func = []( const VectorXcd x ) {
    const VectorXcd res = 2. * ( x.array().pow( 3.0 ) ).matrix() + 10. * x.array().cos().matrix();
    return res;
  };

  VectorXd X( 3 );
  X << 1., 2., 3.;

  auto [F_forward, J_forward] = Complex::forwardDifference( vector_func, X );
  auto J_central              = Complex::centralDifference( vector_func, X );
  auto J_4thorderAccurate     = Complex::fourthOrderAccurateDerivative( vector_func, X );

  auto     df_dx_  = 2. * 3. * X.array().pow( 2.0 ).matrix() + 10. * ( -X.array().sin().matrix() );
  auto     Ftarget = 2. * ( X.array().pow( 3.0 ) ).matrix() + 10. * X.array().cos().matrix();
  MatrixXd Jtarget = MatrixXd::Zero( 3, 3 );
  Jtarget( 0, 0 )  = df_dx_( 0 );
  Jtarget( 1, 1 )  = df_dx_( 1 );
  Jtarget( 2, 2 )  = df_dx_( 2 );

  throwExceptionOnFailure( checkIfEqual< double >( F_forward, Ftarget ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with complex forward difference doesn't yield the right "
                                           "function value" );
  throwExceptionOnFailure( checkIfEqual< double >( J_forward, Jtarget, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with complex forward difference doesn't yield the right "
                                           "derivative" );
  throwExceptionOnFailure( checkIfEqual< double >( J_central, Jtarget, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with complex central difference doesn't yield the right "
                                           "derivative" );
  throwExceptionOnFailure( checkIfEqual< double >( J_4thorderAccurate, Jtarget, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed: jacobian with complex 4th order accuracy doesn't yield the right "
                                           "derivative" );
}

int main()
{
  testNumericDifferentiationForScalars();
  testNumericalDifferntiationForVectorValuedFunctions();
  testComplexStepDerivativeForScalars();
  testComplexStepDerivativeForVectorValuedFunctions();
  return 0;
}
