#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTimeIntegration.h"

using Marmot::MakeString;
using namespace Marmot::NumericalAlgorithms::TimeIntegration;
using namespace Marmot::Testing;

void test_explicitEuler()
{
  // f(x)  = x^2
  // f'(x) = 2x

  auto f     = []( double x ) { return x * x; };
  auto fRate = []( double yN, double x ) { return 2 * x; };

  double x0 = 1; // start = 1; end = 2
  double dx = 1e-6;
  double f0 = f( x0 );

  double f_approx = f0;
  double x        = x0;

  while ( abs( x - 2 ) > 1e-10 ) {
    f_approx = explicitEuler( f_approx, dx, fRate, x );
    x += dx;
  }

  // f_approx should be f(2)
  double computed = f_approx;

  // Expected: f(2) = 2^2 = 4
  double expected = 4.0;

  throwExceptionOnFailure( checkIfEqual( computed, expected, 1e-5 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_semiImplicitEuler()
{
  {
    // Testcase 1
    // f(x)  = x^2
    // f'(x) = 2x

    auto f     = []( double x ) { return x * x; };
    auto fRate = []( const Eigen::VectorXd& yN, double x ) -> Eigen::VectorXd {
      Eigen::VectorXd result( yN.size() );
      result( 0 ) = 2 * x;
      return result;
    };

    double          x0 = 1; // start = 1; end = 2
    double          dx = 1e-6;
    Eigen::VectorXd f0( 1 );
    f0( 0 ) = f( x0 );

    Eigen::VectorXd f_approx = f0;
    double          x        = x0;

    while ( abs( x - 2 ) > 1e-10 ) {
      f_approx = semiImplicitEuler< 1 >( f_approx, dx, fRate, x );
      x += dx;
    }

    // f_approx should be f(2)
    double computed = f_approx( 0 );

    // Expected: f(2) = 2^2 = 4
    double expected = 4.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected, 1e-5 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Testcase 2
    // Define the system of differential equations Dy/dt = A * y
    // where A is a diagonal matrix with entries 1, 2, and 3.
    Eigen::Matrix3d A;
    // clang-format off
    A << 1, 0, 0,
         0, 2, 0,
         0, 0, 3;
    // clang-format on

    auto fRate = [A]( const Eigen::Vector3d& y ) -> Eigen::Vector3d { return A * y; };

    // Initial condition y0 = [1, 1, 1]^T.
    Eigen::Vector3d y0;
    y0 << 1.0, 1.0, 1.0;

    double          dt = 1e-6; // Total time t = steps * dt = 1e6 * 1e-6 = 1.0
    Eigen::Vector3d y  = y0;

    int i = 0;
    while ( i < 1e6 ) {
      y = semiImplicitEuler< 3 >( y, dt, fRate );
      i++;
    }

    // Expected solution after time t = 1.0.
    // For a diagonal matrix A, the solution is y(t) = exp(At) * y0.
    // Since A is diagonal, exp(At) is simply the exponential of each diagonal entry.
    Eigen::Vector3d expected;
    expected( 0 ) = std::exp( 1.0 ); // e^1
    expected( 1 ) = std::exp( 2.0 ); // e^2
    expected( 2 ) = std::exp( 3.0 ); // e^3

    throwExceptionOnFailure( checkIfEqual< double >( y, expected, 1e-4 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }
}

void test_explicitEulerRichardson()
{
  // f(x)  = x^2
  // f'(x) = 2x

  auto f     = []( double x ) { return x * x; };
  auto fRate = []( double yN, double x ) { return 2 * x; };

  double x0 = 1; // start = 1; end = 2
  double dx = 1e-6;
  double f0 = f( x0 );

  double f_approx = f0;
  double x        = x0;

  while ( abs( x - 2 ) > 1e-10 ) {
    f_approx = explicitEulerRichardson( f_approx, dx, fRate, x );
    x += dx;
  }

  // f_approx should be f(2)
  double computed = f_approx;

  // Expected: f(2) = 4
  double expected = 4.0;

  throwExceptionOnFailure( checkIfEqual( computed, expected, 1e-5 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_explicitEulerRichardsonWithErrorEstimator()
{
  // f(x)  = x^2
  // f'(x) = 2x

  auto f     = []( double x ) { return x * x; };
  auto fRate = []( Eigen::Matrix< double, 1, 1 > yN, double x ) {
    Eigen::Matrix< double, 1, 1 > result;
    result( 0 ) = 2 * x;
    return result;
  };

  double x0  = 1.0;   // start = 1.0; end = 2.0
  double dx  = 1e-6;
  double TOL = 1e-12; // tolerance for error estimation
  double f0  = f( x0 );

  Eigen::Matrix< double, 1, 1 > f_approx;
  f_approx( 0 )   = f0; // initial value
  double x        = x0;
  double stepSize = dx;

  while ( abs( x - 2 ) > 1e-15 ) {
    auto [yNew, tauNew] = explicitEulerRichardsonWithErrorEstimator< 1 >( f_approx, stepSize, TOL, fRate, x );
    f_approx( 0 )       = yNew( 0 );

    x += stepSize;
    stepSize = std::min( std::min( stepSize, tauNew ), 2. - x );
  }

  // f_approx should be f(2)
  double computed = f_approx( 0 );

  // Expected: f(2) = 4
  double expected = 4.0;

  throwExceptionOnFailure( checkIfEqual( computed, expected, 1e-6 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ test_explicitEuler,
                                                       test_semiImplicitEuler,
                                                       test_explicitEulerRichardson,
                                                       test_explicitEulerRichardsonWithErrorEstimator };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
