#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Math;
using namespace Marmot::Testing;

void test_linearInterpolation()
{
  // Testcase 1: Normal case
  {
    double x  = 2.0;
    double x0 = 1.0, x1 = 3.0;
    double y0 = 2.0, y1 = 4.0;
    double expected = 3.0; // y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    double computed = linearInterpolation( x, x0, x1, y0, y1 );

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  // Testcase 2: x equals x0
  {
    double x  = 1.0;
    double x0 = 1.0, x1 = 3.0;
    double y0 = 2.0, y1 = 4.0;
    double expected = 2.0; // y0
    double computed = linearInterpolation( x, x0, x1, y0, y1 );

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  // Testcase 3: x equals x1
  {
    double x  = 3.0;
    double x0 = 1.0, x1 = 3.0;
    double y0 = 2.0, y1 = 4.0;
    double expected = 4.0; // y1
    double computed = linearInterpolation( x, x0, x1, y0, y1 );

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }

  // Testcase 4: Negative values
  {
    double x  = -1.0;
    double x0 = -2.0, x1 = 0.0;
    double y0 = -4.0, y1 = 0.0;
    double expected = -2.0;
    double computed = linearInterpolation( x, x0, x1, y0, y1 );

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 4" );
  }

  // Testcase 5: x is outside the range (extrapolation)
  {
    double x  = 5.0;
    double x0 = 1.0, x1 = 3.0;
    double y0 = 2.0, y1 = 4.0;
    double expected = 6.0;
    double computed = linearInterpolation( x, x0, x1, y0, y1 );

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 5" );
  }
}

void test_exp()
{
  // Testcase 1: x = 0, exp(0) = 1
  {
    double computed = Marmot::Math::exp( 0.0 );
    double expected = 1.0;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  // Testcase 2: x = -64, should return 0.0
  {
    double computed = Marmot::Math::exp( -64.0 );
    double expected = 0.0;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  // Testcase 3: x = 64, should return exp(64)
  {
    double computed = Marmot::Math::exp( 64.0 );
    double expected = 6.235149080811617e+27;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }

  // Testcase 4: x = -100, should return 0.0 (underflow)
  {
    double computed = Marmot::Math::exp( -100.0 );
    double expected = 0.0;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 4" );
  }

  // Testcase 5: x = 100, should return exp(64) (overflow)
  {
    double computed = Marmot::Math::exp( 100.0 );
    double expected = 6.235149080811617e+27;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 5" );
  }

  // Testcase 6: x = 10, within bounds
  {
    double computed = Marmot::Math::exp( 10.0 );
    double expected = 22026.465794806718;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 6" );
  }

  // Testcase 7: x = -10, within bounds
  {
    double computed = Marmot::Math::exp( -10.0 );
    double expected = 0.0000453999298;
    throwExceptionOnFailure( checkIfEqual( computed, expected, 1e-13 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 7" );
  }
}

void test_makeReal()
{
  // Testcase 1: Input is a complex number
  {
    std::complex< double > complexInput( 3.0, 4.0 );
    double                 computed = makeReal( complexInput );
    double                 expected = 3.0; // real part of the complex number
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  // Testcase 2: Input is an autodiff::real
  {
    autodiff::real realInput = 5.5;
    double         computed  = makeReal( realInput );
    double         expected  = 5.5;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  // Testcase 3: Input is an autodiff::dual
  {
    autodiff::dual dualInput = 7.7;
    double         computed  = makeReal( dualInput );
    double         expected  = 7.7;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }

  // Testcase 4: Input is a plain double
  {
    double doubleInput = 9.9;
    double computed    = makeReal( doubleInput );
    double expected    = 9.9;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 4" );
  }

  // Testcase 5: Input is an Eigen::Matrix with autodiff::dual elements
  {
    Eigen::Matrix< autodiff::dual, 2, 2 > matrixInput;
    matrixInput << autodiff::dual( 1.1 ), autodiff::dual( 2.2 ), autodiff::dual( 3.3 ), autodiff::dual( 4.4 );

    Eigen::Matrix< double, 2, 2 > computed = makeReal( matrixInput );

    Eigen::Matrix< double, 2, 2 > expected;
    // clang-format off
    expected << 1.1, 2.2,
                3.3, 4.4;
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 5" );
  }

  // Testcase 6: Input is an Eigen::Vector with autodiff::real elements
  {
    Eigen::Vector< autodiff::real, Eigen::Dynamic > vectorInput( 3 );
    vectorInput << autodiff::real( 1.5 ), autodiff::real( 2.5 ), autodiff::real( 3.5 );

    Eigen::VectorXd computed = makeReal( vectorInput );

    Eigen::VectorXd expected( 3 );
    expected << 1.5, 2.5, 3.5;

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 6" );
  }

  // Testcase 7: Input is an Eigen::Matrix with double elements
  {
    Eigen::Matrix< double, 3, 3 > matrixInput;
    // clang-format off
    matrixInput << 1.1, 2.2, 3.3,
                   4.4, 5.5, 6.6,
                   7.7, 8.8, 9.9;
    // clang-format on

    Eigen::Matrix< double, 3, 3 > computed = makeReal( matrixInput );

    Eigen::Matrix< double, 3, 3 > expected;
    // clang-format off
    expected << 1.1, 2.2, 3.3,
                4.4, 5.5, 6.6,
                7.7, 8.8, 9.9;
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 7" );
  }
}

void test_getExponentPowerTen()
{
  // Testcase 1: positive numbers
  {
    double computed = getExponentPowerTen( 5e5 );
    double expected = 5;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }
  {
    double computed = getExponentPowerTen( 1000 );
    double expected = 3;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }
  {
    double computed = getExponentPowerTen( 0.0000001 );
    double expected = -7;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  // Testcase 2: negative numbers
  {
    double computed = getExponentPowerTen( -7e4 );
    double expected = 4;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }
  {
    double computed = getExponentPowerTen( -0.0001 );
    double expected = -4;
    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }
}

void test_orthonormalCoordinateSystem()
{
  {
    // Test 1: Single vector input
    Marmot::Vector3d normalVector( 2.0, 0.0, 1.0 );
    Marmot::Matrix3d computed = orthonormalCoordinateSystem( normalVector );

    Marmot::Matrix3d expected;
    // clang-format off
    expected << 2.0/sqrt(5.0), 0.0, -1.0/sqrt(5.0), // First column: normal vector
                0.0,           1.0, 0.0,            // Second column: orthogonal vector
                1.0/sqrt(5.0), 0.0, 2.0/sqrt(5.0);  // Third column: cross product of the first and second vectors
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Test 2: Two orthogonal vector inputs
    Marmot::Vector3d n1( 1.0, 0.0, 0.0 );
    Marmot::Vector3d n2( 0.0, 1.0, 0.0 );

    Marmot::Matrix3d computed = orthonormalCoordinateSystem( n1, n2 );

    Marmot::Matrix3d expected;
    // clang-format off
    expected << 1.0, 0.0, 0.0,   // First column: n1
                0.0, 1.0, 0.0,   // Second column: n2
                0.0, 0.0, 1.0;   // Third column: cross product of n1 and n2
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }
}

void test_directionCosines()
{
  // Rotation around z-axis
  Marmot::Matrix3d transformedCoordinateSystem;
  // clang-format off
  transformedCoordinateSystem <<  1/sqrt(2), -1/sqrt(2), 0, // First vector
                                  1/sqrt(2),  1/sqrt(2), 0, // Śecond vector
                                  0,          0,         1; // Third vector
  // clang-format on

  Marmot::Matrix3d computed = directionCosines( transformedCoordinateSystem );

  // Expected matrix (rows here correspond to input columns above)
  Marmot::Matrix3d expected;
  // clang-format off
  expected <<  1/sqrt(2), 1/sqrt(2), 0,
              -1/sqrt(2), 1/sqrt(2), 0,
               0,         0,         1;
  // clang-format on

  throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_isNaN()
{
  // Testcase 1: NaN-value
  double nanValue = std::numeric_limits< double >::quiet_NaN();
  bool   computed = isNaN( nanValue );
  bool   expected = true;
  throwExceptionOnFailure( checkIfEqual( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );

  // Testcase 2: Normal value
  double normalValue = 42.0;
  computed           = isNaN( normalValue );
  expected           = false;
  throwExceptionOnFailure( checkIfEqual( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );

  // Testcase 3: Infinity
  double infValue = std::numeric_limits< double >::infinity();
  computed        = isNaN( infValue );
  expected        = false;
  throwExceptionOnFailure( checkIfEqual( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );

  // Testcase 4: -Infinity
  double negInfValue = -std::numeric_limits< double >::infinity();
  computed           = isNaN( negInfValue );
  expected           = false;
  throwExceptionOnFailure( checkIfEqual( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 4" );

  // Testcase 5: zero value
  double zeroValue = 0.0;
  computed         = isNaN( zeroValue );
  expected         = false;
  throwExceptionOnFailure( checkIfEqual( computed, expected ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 5" );
}

void test_radToDeg()
{
  double alpha    = Marmot::Constants::Pi / 4; // 45 degrees in radians
  double computed = radToDeg( alpha );
  double expected = 45.0;

  throwExceptionOnFailure( checkIfEqual( computed, expected ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_degToRad()
{
  double alpha    = 90.0;                      // 90 degrees
  double computed = degToRad( alpha );
  double expected = Marmot::Constants::Pi / 2; // radians

  throwExceptionOnFailure( checkIfEqual( computed, expected ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_macauly()
{
  {
    // Testcase 1: Positive value
    double positiveValue = 5.0;
    double computed      = macauly( positiveValue );
    double expected      = 5.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Testcase 2: Zero value
    double zeroValue = 0.0;
    double computed  = macauly( zeroValue );
    double expected  = 0.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  {
    // Testcase 3: Negative value
    double negativeValue = -3.0;
    double computed      = macauly( negativeValue );
    double expected      = 0.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }
}

void test_macaulyMatrix()
{
  {
    // Testcase 1: Matrix with positive and negative values
    Eigen::Matrix2d inputMatrix;
    // clang-format off
    inputMatrix << 1.0, -2.0,
                  -3.0,  4.0;
    // clang-format on

    Eigen::Matrix2d computed = macaulyMatrix( inputMatrix );

    Eigen::Matrix2d expected;
    // clang-format off
    expected << 1.0, 0.0,
                0.0, 4.0;
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Testcase 2: Zero matrix
    Eigen::Matrix2d inputMatrix = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d computed    = macaulyMatrix( inputMatrix );

    Eigen::Matrix2d expected = Eigen::Matrix2d::Zero();

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  {
    // Testcase 3: Matrix with all positive values
    Eigen::Matrix2d inputMatrix;
    // clang-format off
    inputMatrix << 5.0, 6.0,
                   7.0, 8.0;
    // clang-format on

    Eigen::Matrix2d computed = macaulyMatrix( inputMatrix );

    Eigen::Matrix2d expected;
    // clang-format off
    expected << 5.0, 6.0,
                7.0, 8.0;
    // clang-format on

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }

  {
    // Testcase 4: Matrix with all negative values
    Eigen::Matrix2d inputMatrix;
    // clang-format off
    inputMatrix << -1.0, -2.0,
                   -3.0, -4.0;
    // clang-format on

    Eigen::Matrix2d computed = macaulyMatrix( inputMatrix );

    Eigen::Matrix2d expected = Eigen::Matrix2d::Zero();

    throwExceptionOnFailure( checkIfEqual< double >( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 4" );
  }
}

void test_heaviside()
{
  {
    // Testcase 1: Positive value
    double positiveValue = 5.0;
    double computed      = heaviside( positiveValue );
    int    expected      = 1;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Testcase 2: Zero value
    double zeroValue = 0.0;
    double computed  = heaviside( zeroValue );
    int    expected  = 1;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  {
    // Testcase 3: Negative value
    double negativeValue = -3.0;
    double computed      = heaviside( negativeValue );
    int    expected      = 0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }
}

void test_sgn()
{
  {
    // Testcase 1: Positive value
    double positiveValue = 5.0;
    double computed      = sgn( positiveValue );
    double expected      = 1.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  {
    // Testcase 2: Zero value
    double zeroValue = 0.0;
    double computed  = sgn( zeroValue );
    double expected  = 0.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  {
    // Testcase 3: Negative value
    double negativeValue = -3.0;
    double computed      = sgn( negativeValue );
    double expected      = -1.0;

    throwExceptionOnFailure( checkIfEqual( computed, expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }
}

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
  auto tests = std::vector< std::function< void() > >{ test_linearInterpolation,
                                                       test_exp,
                                                       test_makeReal,
                                                       test_getExponentPowerTen,
                                                       test_orthonormalCoordinateSystem,
                                                       test_directionCosines,
                                                       test_isNaN,
                                                       test_radToDeg,
                                                       test_degToRad,
                                                       test_macauly,
                                                       test_macaulyMatrix,
                                                       test_heaviside,
                                                       test_sgn,
                                                       test_explicitEuler,
                                                       test_semiImplicitEuler,
                                                       test_explicitEulerRichardson,
                                                       test_explicitEulerRichardsonWithErrorEstimator };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
