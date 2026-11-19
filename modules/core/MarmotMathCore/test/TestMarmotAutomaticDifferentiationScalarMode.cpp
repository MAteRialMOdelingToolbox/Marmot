#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <cmath>
#include <functional>
#include <numbers>

using namespace Marmot::Testing;
using namespace Marmot::Math::Differentiation::Automatic::ScalarMode;

namespace {

  void check( double a, double b, const std::string& what, double tol = 1e-12 )
  {
    throwExceptionOnFailure( checkIfEqual( a, b, tol * std::max( 1.0, std::abs( b ) ) ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " " << what << ": got " << a << ", expected " << b );
  }

  /** f, f', f'', f''' of a function at x, for one test case */
  struct Case {
    std::string                                       name;
    std::function< HyperDual< 3 >( HyperDual< 3 > ) > f3;
    std::function< Taylor< 3 >( Taylor< 3 > ) >       taylor;
    std::array< double, 4 >                           derivatives;
  };

} // namespace

void testHyperDualArithmetic()
{
  // f( x, y ) = x^2 y / ( x + y ): mixed derivatives from one hyper-dual evaluation
  using D2       = HyperDual< 2 >;
  const double x = 0.7, y = 1.3;
  D2           x_ = x, y_ = y;
  x_.grad()       = 1.0; // level 1: x
  y_.val().grad() = 1.0; // level 2: y
  const D2     f  = x_ * x_ * y_ / ( x_ + y_ );
  const double s  = x + y;
  check( double( f ), x * x * y / s, "value" );
  check( primalValue( f.grad() ), ( 2 * x * y * s - x * x * y ) / ( s * s ), "df/dx" );
  check( f.val().grad(), ( x * x * s - x * x * y ) / ( s * s ), "df/dy" );
  // d2f/dxdy = d/dx [ x^3 / s^2 ]
  const double d2 = x * x * ( x + 3 * y ) / ( s * s * s );
  check( f.grad().grad(), d2, "d2f/dxdy" );

  // operators with constants and compound assignment
  D2 g = x_;
  g += 2.;
  g *= 3;
  g -= x_;
  g /= 2.;
  check( double( g ), ( ( x + 2 ) * 3 - x ) / 2, "compound value" );
  check( primalValue( g.grad() ), 1.0, "compound derivative" );
  check( double( 1. - x_ ), 1 - x, "constant minus" );
  check( primalValue( ( 2. / x_ ).grad() ), -2 / ( x * x ), "constant over" );
  throwExceptionOnFailure( x_ < y_ && x_ < 1 && 0 < x_ && x_ != y_ && x_ == x && max( x_, y_ ) == y_ &&
                             min( x_, 0.5 ) == 0.5,
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " comparisons failed" );
}

void testElementaryFunctions()
{
  const double x = 0.6;
  // f, f', f'', f''' at x
  const std::vector< Case > cases = {
    { "sqrt",
      []( HyperDual< 3 > a ) { return sqrt( a ); },
      []( Taylor< 3 > a ) { return sqrt( a ); },
      { std::sqrt( x ), 0.5 / std::sqrt( x ), -0.25 * std::pow( x, -1.5 ), 0.375 * std::pow( x, -2.5 ) } },
    { "cbrt",
      []( HyperDual< 3 > a ) { return cbrt( a ); },
      []( Taylor< 3 > a ) { return cbrt( a ); },
      { std::cbrt( x ),
        1. / 3 * std::pow( x, -2. / 3 ),
        -2. / 9 * std::pow( x, -5. / 3 ),
        10. / 27 * std::pow( x, -8. / 3 ) } },
    { "exp",
      []( HyperDual< 3 > a ) { return exp( a ); },
      []( Taylor< 3 > a ) { return exp( a ); },
      { std::exp( x ), std::exp( x ), std::exp( x ), std::exp( x ) } },
    { "log",
      []( HyperDual< 3 > a ) { return log( a ); },
      []( Taylor< 3 > a ) { return log( a ); },
      { std::log( x ), 1 / x, -1 / ( x * x ), 2 / ( x * x * x ) } },
    { "log10",
      []( HyperDual< 3 > a ) { return log10( a ); },
      []( Taylor< 3 > a ) { return log10( a ); },
      { std::log10( x ),
        1 / ( x * std::numbers::ln10 ),
        -1 / ( x * x * std::numbers::ln10 ),
        2 / ( x * x * x * std::numbers::ln10 ) } },
    { "pow( x, 2.5 )",
      []( HyperDual< 3 > a ) { return pow( a, 2.5 ); },
      []( Taylor< 3 > a ) { return pow( a, 2.5 ); },
      { std::pow( x, 2.5 ), 2.5 * std::pow( x, 1.5 ), 3.75 * std::pow( x, 0.5 ), 1.875 * std::pow( x, -0.5 ) } },
    { "pow( x, 2 )",
      []( HyperDual< 3 > a ) { return pow( a, 2 ); },
      []( Taylor< 3 > a ) { return pow( a, 2 ); },
      { x * x, 2 * x, 2, 0 } },
    { "pow( 2, x )",
      []( HyperDual< 3 > a ) { return pow( 2., a ); },
      []( Taylor< 3 > a ) { return pow( 2., a ); },
      { std::pow( 2, x ),
        std::pow( 2, x ) * std::log( 2 ),
        std::pow( 2, x ) * std::pow( std::log( 2 ), 2 ),
        std::pow( 2, x ) * std::pow( std::log( 2 ), 3 ) } },
    { "pow( x, x )",
      []( HyperDual< 3 > a ) { return pow( a, a ); },
      []( Taylor< 3 > a ) { return pow( a, a ); },
      { std::pow( x, x ),
        std::pow( x, x ) * ( std::log( x ) + 1 ),
        std::pow( x, x ) * ( std::pow( std::log( x ) + 1, 2 ) + 1 / x ),
        std::pow( x, x ) * ( std::pow( std::log( x ) + 1, 3 ) + 3 * ( std::log( x ) + 1 ) / x - 1 / ( x * x ) ) } },
    { "sin",
      []( HyperDual< 3 > a ) { return sin( a ); },
      []( Taylor< 3 > a ) { return sin( a ); },
      { std::sin( x ), std::cos( x ), -std::sin( x ), -std::cos( x ) } },
    { "cos",
      []( HyperDual< 3 > a ) { return cos( a ); },
      []( Taylor< 3 > a ) { return cos( a ); },
      { std::cos( x ), -std::sin( x ), -std::cos( x ), std::sin( x ) } },
    { "tan",
      []( HyperDual< 3 > a ) { return tan( a ); },
      []( Taylor< 3 > a ) { return tan( a ); },
      { std::tan( x ),
        1 + std::pow( std::tan( x ), 2 ),
        2 * std::tan( x ) * ( 1 + std::pow( std::tan( x ), 2 ) ),
        2 * ( 1 + std::pow( std::tan( x ), 2 ) ) * ( 1 + 3 * std::pow( std::tan( x ), 2 ) ) } },
    { "asin",
      []( HyperDual< 3 > a ) { return asin( a ); },
      []( Taylor< 3 > a ) { return asin( a ); },
      { std::asin( x ),
        std::pow( 1 - x * x, -0.5 ),
        x * std::pow( 1 - x * x, -1.5 ),
        ( 1 + 2 * x * x ) * std::pow( 1 - x * x, -2.5 ) } },
    { "acos",
      []( HyperDual< 3 > a ) { return acos( a ); },
      []( Taylor< 3 > a ) { return acos( a ); },
      { std::acos( x ),
        -std::pow( 1 - x * x, -0.5 ),
        -x * std::pow( 1 - x * x, -1.5 ),
        -( 1 + 2 * x * x ) * std::pow( 1 - x * x, -2.5 ) } },
    { "atan",
      []( HyperDual< 3 > a ) { return atan( a ); },
      []( Taylor< 3 > a ) { return atan( a ); },
      { std::atan( x ),
        1 / ( 1 + x * x ),
        -2 * x / std::pow( 1 + x * x, 2 ),
        ( 6 * x * x - 2 ) / std::pow( 1 + x * x, 3 ) } },
    { "sinh",
      []( HyperDual< 3 > a ) { return sinh( a ); },
      []( Taylor< 3 > a ) { return sinh( a ); },
      { std::sinh( x ), std::cosh( x ), std::sinh( x ), std::cosh( x ) } },
    { "cosh",
      []( HyperDual< 3 > a ) { return cosh( a ); },
      []( Taylor< 3 > a ) { return cosh( a ); },
      { std::cosh( x ), std::sinh( x ), std::cosh( x ), std::sinh( x ) } },
    { "tanh",
      []( HyperDual< 3 > a ) { return tanh( a ); },
      []( Taylor< 3 > a ) { return tanh( a ); },
      { std::tanh( x ),
        1 - std::pow( std::tanh( x ), 2 ),
        -2 * std::tanh( x ) * ( 1 - std::pow( std::tanh( x ), 2 ) ),
        ( 1 - std::pow( std::tanh( x ), 2 ) ) * ( 6 * std::pow( std::tanh( x ), 2 ) - 2 ) } },
    { "abs( -x )",
      []( HyperDual< 3 > a ) { return abs( -a ); },
      []( Taylor< 3 > a ) { return abs( -a ); },
      { x, 1, 0, 0 } },
    { "1 / ( 1 + x^2 ) - x",
      []( HyperDual< 3 > a ) { return 1. / ( 1. + a * a ) - a; },
      []( Taylor< 3 > a ) { return 1. / ( 1. + a * a ) - a; },
      { 1 / ( 1 + x * x ) - x,
        -2 * x / std::pow( 1 + x * x, 2 ) - 1,
        ( 6 * x * x - 2 ) / std::pow( 1 + x * x, 3 ),
        24 * x * ( 1 - x * x ) / std::pow( 1 + x * x, 4 ) } },
  };

  for ( const auto& c : cases ) {
    // the same variable on all three levels of a 3rd order hyper-dual number
    HyperDual< 3 > x3       = x;
    x3.grad()               = 1.0;
    x3.val().grad()         = 1.0;
    x3.val().val().grad()   = 1.0;
    const HyperDual< 3 > y3 = c.f3( x3 );
    check( double( y3 ), c.derivatives[0], c.name + " value (hyper-dual)" );
    check( primalValue( y3.grad() ), c.derivatives[1], c.name + " 1st derivative (hyper-dual)" );
    check( primalValue( y3.grad().grad() ), c.derivatives[2], c.name + " 2nd derivative (hyper-dual)" );
    check( y3.grad().grad().grad(), c.derivatives[3], c.name + " 3rd derivative (hyper-dual)" );

    Taylor< 3 > t               = x;
    t[1]                        = 1.0;
    const Taylor< 3 > yt        = c.taylor( t );
    double            factorial = 1.0;
    for ( std::size_t k = 0; k <= 3; k++ ) {
      factorial *= k ? double( k ) : 1.0;
      check( yt[k] * factorial, c.derivatives[k], c.name + " derivative " + std::to_string( k ) + " (Taylor)" );
    }
  }
}

void testTaylorHighOrder()
{
  // exp( -2 t ) and t^-0.3 up to the 7th derivative, as the Post-Widder formulas of the creep models need them
  const double t              = 1.7;
  Taylor< 7 >  t_             = t;
  t_[1]                       = 1.0;
  const Taylor< 7 > e         = exp( -2. * t_ );
  const Taylor< 7 > p         = pow( t_, -0.3 );
  double            factorial = 1.0, binomial = 1.0;
  for ( std::size_t k = 0; k <= 7; k++ ) {
    factorial *= k ? double( k ) : 1.0;
    check( e[k] * factorial,
           std::pow( -2.0, double( k ) ) * std::exp( -2 * t ),
           "exp( -2 t ), order " + std::to_string( k ) );
    check( p[k], binomial * std::pow( t, -0.3 - double( k ) ), "t^-0.3, order " + std::to_string( k ) );
    binomial *= ( -0.3 - double( k ) ) / double( k + 1 );
  }
}

void testTaylorAtZero()
{
  // a function of t that is identically zero, e.g. B4's drying compliance for a vanishing humidity difference: sqrt and
  // pow of it must be zero, not NaN
  Taylor< 5 > t           = 0.3;
  t[1]                    = 1.0;
  const Taylor< 5 > zero  = exp( t * 0.0 ) - 1.0;
  const Taylor< 5 > root  = sqrt( zero );
  const Taylor< 5 > power = pow( zero, 1.5 );
  for ( std::size_t k = 0; k <= 5; k++ ) {
    check( root[k], 0.0, "sqrt( 0 ), order " + std::to_string( k ) );
    check( power[k], 0.0, "pow( 0, 1.5 ), order " + std::to_string( k ) );
  }
}

void testEigen()
{
  // gradient of |A x|^2 / 2 with dual numbers in Eigen, one direction per evaluation
  using D = HyperDual< 1 >;
  Marmot::Math::Matrix3d A;
  A << 2, 1, 0, 1, 3, 1, 0, 1, 4;
  Marmot::Math::Vector3d       x( 0.1, -0.2, 0.3 );
  const Marmot::Math::Vector3d g = A.transpose() * A * x;
  for ( int j = 0; j < 3; j++ ) {
    Eigen::Vector< D, 3 > x_      = x.cast< D >();
    x_( j ).grad()                = 1.0;
    const Eigen::Vector< D, 3 > y = A * x_;
    const D                     f = 0.5 * y.squaredNorm();
    check( f.grad(), g( j ), "Eigen gradient" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testHyperDualArithmetic,
                                                       testElementaryFunctions,
                                                       testTaylorHighOrder,
                                                       testTaylorAtZero,
                                                       testEigen };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
