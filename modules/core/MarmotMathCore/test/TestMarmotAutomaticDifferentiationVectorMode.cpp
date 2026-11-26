#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorModeSecondOrder.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorModeThirdOrder.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <cmath>
#include <functional>
#include <numbers>

using namespace Marmot::Testing;
using Marmot::Math::Differentiation::Automatic::VectorMode::Dual;

namespace {

  using Dual2 = Dual< 2 >;

  /** x and y as the independent variables 0 and 1 */
  std::pair< Dual2, Dual2 > independentVariables( double x, double y )
  {
    Dual2 x_( x ), y_( y );
    x_.dx( 0 ) = 1.0;
    y_.dx( 1 ) = 1.0;
    return { x_, y_ };
  }

  void checkDual( const Dual2&       f,
                  double             value,
                  double             df_dx,
                  double             df_dy,
                  const std::string& what,
                  double             tol = 1e-13 )
  {
    throwExceptionOnFailure( checkIfEqual( f.val(), value, tol ) && checkIfEqual( f.dx( 0 ), df_dx, tol ) &&
                               checkIfEqual( f.dx( 1 ), df_dy, tol ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " " << what << ": got " << f << ", expected " << value << " ["
                               << df_dx << ", " << df_dy << "]" );
  }

} // namespace

void testArithmetic()
{
  const double x = 0.7, y = -1.3;
  auto [x_, y_] = independentVariables( x, y );

  checkDual( x_ + y_, x + y, 1, 1, "x + y" );
  checkDual( x_ - y_, x - y, 1, -1, "x - y" );
  checkDual( x_ * y_, x * y, y, x, "x * y" );
  checkDual( x_ / y_, x / y, 1 / y, -x / ( y * y ), "x / y" );
  checkDual( -x_, -x, -1, 0, "-x" );
  checkDual( +x_, x, 1, 0, "+x" );
  checkDual( x_ + 2, x + 2, 1, 0, "x + 2" );
  checkDual( 2. + x_, x + 2, 1, 0, "2 + x" );
  checkDual( x_ - 2, x - 2, 1, 0, "x - 2" );
  checkDual( 2. - x_, 2 - x, -1, 0, "2 - x" );
  checkDual( 3 * x_, 3 * x, 3, 0, "3 * x" );
  checkDual( x_ * 3., 3 * x, 3, 0, "x * 3" );
  checkDual( x_ / 4., x / 4, 0.25, 0, "x / 4" );
  checkDual( 4. / x_, 4 / x, -4 / ( x * x ), 0, "4 / x" );

  // nested expressions, with intermediate nodes kept by value
  checkDual( ( x_ * y_ + 2. * x_ ) / ( y_ - x_ * x_ ),
             ( x * y + 2 * x ) / ( y - x * x ),
             ( ( y + 2 ) * ( y - x * x ) + ( x * y + 2 * x ) * 2 * x ) / ( ( y - x * x ) * ( y - x * x ) ),
             ( x * ( y - x * x ) - ( x * y + 2 * x ) ) / ( ( y - x * x ) * ( y - x * x ) ),
             "( x y + 2 x ) / ( y - x^2 )" );
}

void testAssignmentAndAliasing()
{
  const double x = 0.7, y = -1.3;
  auto [x_, y_] = independentVariables( x, y );

  // the right-hand side refers to the left-hand side
  Dual2 z = x_;
  z       = z * y_ + z;
  checkDual( z, x * y + x, y + 1, x, "z = z * y + z" );

  z = x_;
  z += y_;
  checkDual( z, x + y, 1, 1, "z += y" );
  z -= 2. * y_;
  checkDual( z, x - y, 1, -1, "z -= 2 y" );
  z *= x_;
  checkDual( z, ( x - y ) * x, 2 * x - y, -x, "z *= x" );
  z /= y_;
  checkDual( z, ( x - y ) * x / y, ( 2 * x - y ) / y, -x / y - ( x - y ) * x / ( y * y ), "z /= y" );

  z = x_;
  z += 1;
  z -= 0.5;
  z *= 2;
  z /= 4.;
  checkDual( z, ( x + 0.5 ) / 2, 0.5, 0, "constant compound assignment" );

  z = 3.;
  checkDual( z, 3, 0, 0, "assignment of a constant" );

  const Dual2 zero;
  checkDual( zero, 0, 0, 0, "default construction" );

  throwExceptionOnFailure( ( x_ < y_ ) == ( x < y ) && x_ > 0 && 0 < x_ && x_ != y_ && x_ == x && x_ * 1. >= x_ &&
                             x_ <= x_ + y_ * y_,
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " comparisons failed" );
}

void testElementaryFunctions()
{
  const double x = 0.7, y = 1.3;
  auto [x_, y_] = independentVariables( x, y );

  checkDual( sqrt( x_ ), std::sqrt( x ), 0.5 / std::sqrt( x ), 0, "sqrt" );
  checkDual( cbrt( x_ ), std::cbrt( x ), 1. / ( 3 * std::cbrt( x * x ) ), 0, "cbrt" );
  checkDual( exp( x_ ), std::exp( x ), std::exp( x ), 0, "exp" );
  checkDual( log( x_ ), std::log( x ), 1 / x, 0, "log" );
  checkDual( log10( x_ ), std::log10( x ), 1 / ( x * std::numbers::ln10 ), 0, "log10" );
  checkDual( pow( x_, 3 ), std::pow( x, 3 ), 3 * x * x, 0, "pow( x, 3 )" );
  checkDual( pow( x_, -0.5 ), std::pow( x, -0.5 ), -0.5 * std::pow( x, -1.5 ), 0, "pow( x, -0.5 )" );
  checkDual( pow( 2., x_ ), std::pow( 2, x ), std::pow( 2, x ) * std::log( 2 ), 0, "pow( 2, x )" );
  checkDual( pow( x_, y_ ),
             std::pow( x, y ),
             y * std::pow( x, y - 1 ),
             std::pow( x, y ) * std::log( x ),
             "pow( x, y )" );
  checkDual( sin( x_ ), std::sin( x ), std::cos( x ), 0, "sin" );
  checkDual( cos( x_ ), std::cos( x ), -std::sin( x ), 0, "cos" );
  checkDual( tan( x_ ), std::tan( x ), 1 / ( std::cos( x ) * std::cos( x ) ), 0, "tan" );
  checkDual( asin( x_ ), std::asin( x ), 1 / std::sqrt( 1 - x * x ), 0, "asin" );
  checkDual( acos( x_ ), std::acos( x ), -1 / std::sqrt( 1 - x * x ), 0, "acos" );
  checkDual( atan( x_ ), std::atan( x ), 1 / ( 1 + x * x ), 0, "atan" );
  checkDual( atan2( y_, x_ ), std::atan2( y, x ), -y / ( x * x + y * y ), x / ( x * x + y * y ), "atan2" );
  checkDual( sinh( x_ ), std::sinh( x ), std::cosh( x ), 0, "sinh" );
  checkDual( cosh( x_ ), std::cosh( x ), std::sinh( x ), 0, "cosh" );
  checkDual( tanh( x_ ), std::tanh( x ), 1 / ( std::cosh( x ) * std::cosh( x ) ), 0, "tanh" );
  checkDual( abs( -x_ ), x, 1, 0, "abs" );
  checkDual( fabs( x_ - y_ ), y - x, -1, 1, "fabs" );
  checkDual( max( x_, y_ ), y, 0, 1, "max" );
  checkDual( min( x_, y_ ), x, 1, 0, "min" );

  // the chain rule through nested elementary functions
  checkDual( exp( sqrt( x_ * y_ ) ),
             std::exp( std::sqrt( x * y ) ),
             std::exp( std::sqrt( x * y ) ) * 0.5 / std::sqrt( x * y ) * y,
             std::exp( std::sqrt( x * y ) ) * 0.5 / std::sqrt( x * y ) * x,
             "exp( sqrt( x y ) )" );
}

void testEigen()
{
  using Dual6 = Dual< 6 >;

  Marmot::Math::Matrix6d A;
  for ( int i = 0; i < 6; i++ )
    for ( int j = 0; j < 6; j++ )
      A( i, j ) = 1.0 / ( 1 + i + 2 * j );
  Marmot::Math::Vector6d x;
  x << 0.1, -0.2, 0.3, 0.4, -0.5, 0.6;

  Eigen::Vector< Dual6, 6 > x_ = x.cast< Dual6 >();
  for ( int i = 0; i < 6; i++ )
    x_( i ).dx( i ) = 1.0;

  // y = A x, and |y|
  const Eigen::Vector< Dual6, 6 > y_    = A * x_;
  const Dual6                     norm_ = y_.norm();
  const Marmot::Math::Vector6d    y     = A * x;
  const Marmot::Math::Vector6d    dNorm = A.transpose() * y / y.norm();

  for ( int i = 0; i < 6; i++ ) {
    throwExceptionOnFailure( checkIfEqual( y_( i ).val(), y( i ), 1e-15 ),
                             Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " A x" );
    throwExceptionOnFailure( checkIfEqual( norm_.dx( i ), dNorm( i ), 1e-14 ),
                             Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " d|A x|/dx" );
    for ( int j = 0; j < 6; j++ )
      throwExceptionOnFailure( checkIfEqual( y_( i ).dx( j ), A( i, j ), 1e-15 ),
                               Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " d(A x)/dx" );
  }
}

void testFastor()
{
  using Dual9 = Dual< 9 >;
  using namespace Marmot::Math::Tensors::FastorIndices;

  Fastor::Tensor< double, 3, 3 > F;
  F = { { 1.1, 0.2, -0.1 }, { 0.05, 0.9, 0.3 }, { -0.2, 0.1, 1.2 } };

  Fastor::Tensor< Dual9, 3, 3 > F_;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      F_( i, j )                 = F( i, j );
      F_( i, j ).dx( 3 * i + j ) = 1.0;
    }

  // d det( F ) / dF = det( F ) F^-T
  const Dual9                          J_    = Fastor::determinant( F_ );
  const Fastor::Tensor< double, 3, 3 > dJ_dF = Fastor::determinant( F ) * Fastor::transpose( Fastor::inverse( F ) );

  // d tr( F^T F ) / dF = 2 F
  const Fastor::Tensor< Dual9, 3, 3 > C_  = Fastor::einsum< kI, kJ >( F_, F_ );
  const Dual9                         trC = Fastor::trace( C_ );

  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      throwExceptionOnFailure( checkIfEqual( J_.dx( 3 * i + j ), dJ_dF( i, j ), 1e-14 ),
                               Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " d det( F ) / dF" );
      throwExceptionOnFailure( checkIfEqual( trC.dx( 3 * i + j ), 2 * F( i, j ), 1e-14 ),
                               Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " d tr( F^T F ) / dF" );
    }
}

namespace {

  using H3 = Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< 3 >;

  H3 variable( double value, std::size_t direction )
  {
    H3 x = value;
    x.seed( direction );
    return x;
  }

  void checkClose( double a, double b, const std::string& what, double tol = 1e-12 )
  {
    throwExceptionOnFailure( std::abs( a - b ) <= tol * std::max( 1.0, std::abs( b ) ),
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " " << what << ": got " << a << ", expected " << b );
  }

  /** a function of three variables, for second- and third-order vector mode and nested hyper-dual numbers alike */
  template < typename T >
  T composite( const T& x, const T& y, const T& z )
  {
    return exp( x * y ) * sin( z ) / ( 1. + x ) + pow( x, 2.5 ) * log( y ) - sqrt( x * x + z * z ) + atan2( y, z ) -
           tanh( 2. * x - y ) * cbrt( z ) + pow( y, z ) - 3. / ( x + y + z );
  }

} // namespace

void testSecondOrderElementaryFunctions()
{
  // f( x ) with x the variable 1 of 3: value, f', f'' in ( 1, 1 ), everything else zero
  const double x = 0.6;
  struct Case {
    std::string                      name;
    std::function< H3( const H3& ) > f;
    double                           f0, f1, f2;
  };
  const std::vector< Case > cases = {
    { "sqrt",
      []( const H3& a ) { return sqrt( a ); },
      std::sqrt( x ),
      0.5 / std::sqrt( x ),
      -0.25 * std::pow( x, -1.5 ) },
    { "cbrt",
      []( const H3& a ) { return cbrt( a ); },
      std::cbrt( x ),
      std::pow( x, -2. / 3 ) / 3,
      -2. / 9 * std::pow( x, -5. / 3 ) },
    { "exp", []( const H3& a ) { return exp( a ); }, std::exp( x ), std::exp( x ), std::exp( x ) },
    { "log", []( const H3& a ) { return log( a ); }, std::log( x ), 1 / x, -1 / ( x * x ) },
    { "log10",
      []( const H3& a ) { return log10( a ); },
      std::log10( x ),
      1 / ( x * std::numbers::ln10 ),
      -1 / ( x * x * std::numbers::ln10 ) },
    { "pow( x, 2.5 )",
      []( const H3& a ) { return pow( a, 2.5 ); },
      std::pow( x, 2.5 ),
      2.5 * std::pow( x, 1.5 ),
      3.75 * std::pow( x, 0.5 ) },
    { "pow( x, 2 )", []( const H3& a ) { return pow( a, 2 ); }, x * x, 2 * x, 2 },
    { "pow( x, 0 )", []( const H3& a ) { return pow( a, 0 ); }, 1, 0, 0 },
    { "pow( 2, x )",
      []( const H3& a ) { return pow( 2., a ); },
      std::pow( 2, x ),
      std::pow( 2, x ) * std::log( 2 ),
      std::pow( 2, x ) * std::pow( std::log( 2 ), 2 ) },
    { "sin", []( const H3& a ) { return sin( a ); }, std::sin( x ), std::cos( x ), -std::sin( x ) },
    { "cos", []( const H3& a ) { return cos( a ); }, std::cos( x ), -std::sin( x ), -std::cos( x ) },
    { "tan",
      []( const H3& a ) { return tan( a ); },
      std::tan( x ),
      1 + std::pow( std::tan( x ), 2 ),
      2 * std::tan( x ) * ( 1 + std::pow( std::tan( x ), 2 ) ) },
    { "asin",
      []( const H3& a ) { return asin( a ); },
      std::asin( x ),
      std::pow( 1 - x * x, -0.5 ),
      x * std::pow( 1 - x * x, -1.5 ) },
    { "acos",
      []( const H3& a ) { return acos( a ); },
      std::acos( x ),
      -std::pow( 1 - x * x, -0.5 ),
      -x * std::pow( 1 - x * x, -1.5 ) },
    { "atan",
      []( const H3& a ) { return atan( a ); },
      std::atan( x ),
      1 / ( 1 + x * x ),
      -2 * x / std::pow( 1 + x * x, 2 ) },
    { "sinh", []( const H3& a ) { return sinh( a ); }, std::sinh( x ), std::cosh( x ), std::sinh( x ) },
    { "cosh", []( const H3& a ) { return cosh( a ); }, std::cosh( x ), std::sinh( x ), std::cosh( x ) },
    { "tanh",
      []( const H3& a ) { return tanh( a ); },
      std::tanh( x ),
      1 - std::pow( std::tanh( x ), 2 ),
      -2 * std::tanh( x ) * ( 1 - std::pow( std::tanh( x ), 2 ) ) },
    { "abs( -x )", []( const H3& a ) { return abs( -a ); }, x, 1, 0 },
    { "1 / x", []( const H3& a ) { return 1. / a; }, 1 / x, -1 / ( x * x ), 2 / ( x * x * x ) },
    { "2 - 3 x", []( const H3& a ) { return 2. - a * 3.; }, 2 - 3 * x, -3, 0 },
  };
  for ( const auto& c : cases ) {
    const H3 y = c.f( variable( x, 1 ) );
    checkClose( y.val(), c.f0, c.name + " value" );
    checkClose( y.dx( 1 ), c.f1, c.name + " first derivative" );
    checkClose( y.dxx( 1, 1 ), c.f2, c.name + " second derivative" );
    for ( std::size_t i = 0; i < 3; i++ )
      for ( std::size_t j = 0; j < 3; j++ )
        if ( i != 1 || j != 1 )
          checkClose( y.dxx( i, j ), 0.0, c.name + " other second derivatives" );
    throwExceptionOnFailure( y.dx( 0 ) == 0.0 && y.dx( 2 ) == 0.0,
                             Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << c.name << " other derivatives" );
  }
}

void testSecondOrderBinaryOperations()
{
  // x and y are the variables 0 and 2 of 3, so the operands depend on different directions
  const double x = 0.7, y = 1.3;
  const H3     x_ = variable( x, 0 ), y_ = variable( y, 2 );

  const H3 p = x_ * y_;
  checkClose( p.val(), x * y, "x y" );
  checkClose( p.dx( 0 ), y, "d( x y )/dx" );
  checkClose( p.dx( 2 ), x, "d( x y )/dy" );
  checkClose( p.dxx( 0, 2 ), 1, "d2( x y )/dxdy" );
  checkClose( p.dxx( 2, 0 ), 1, "d2( x y )/dydx" );
  checkClose( p.dxx( 0, 0 ) + p.dxx( 2, 2 ) + p.dxx( 1, 1 ), 0, "d2( x y )/dx2, dy2" );

  const H3 q = x_ / y_;
  checkClose( q.dx( 2 ), -x / ( y * y ), "d( x / y )/dy" );
  checkClose( q.dxx( 0, 2 ), -1 / ( y * y ), "d2( x / y )/dxdy" );
  checkClose( q.dxx( 2, 2 ), 2 * x / ( y * y * y ), "d2( x / y )/dy2" );

  const H3 w = pow( x_, y_ );
  checkClose( w.dxx( 0, 2 ), std::pow( x, y - 1 ) * ( 1 + y * std::log( x ) ), "d2 x^y/dxdy" );

  const H3     a  = atan2( y_, x_ );
  const double r2 = x * x + y * y;
  checkClose( a.val(), std::atan2( y, x ), "atan2" );
  checkClose( a.dx( 0 ), -y / r2, "d atan2/dx" );
  checkClose( a.dxx( 0, 2 ), ( y * y - x * x ) / ( r2 * r2 ), "d2 atan2/dxdy" );

  H3 s = x_;
  s += y_;
  s -= 2. * x_;
  s *= y_;
  s /= x_;
  // ( y - x ) y / x
  checkClose( s.val(), ( y - x ) * y / x, "compound assignment" );
  checkClose( s.dxx( 0, 2 ), -2 * y / ( x * x ), "compound assignment d2/dxdy" );
  throwExceptionOnFailure( x_ < y_ && x_ == x && 0 < x_ && s != x_,
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " comparisons" );
}

void testSecondOrderSpecialCases()
{
  // a vanishing base: pow without dividing by the base
  const H3 p = pow( variable( 0.0, 1 ), 2.5 );
  checkClose( p.val(), 0.0, "0^2.5" );
  checkClose( p.dx( 1 ), 0.0, "d x^2.5/dx at 0" );
  checkClose( p.dxx( 1, 1 ), 0.0, "d2 x^2.5/dx2 at 0" );

  // the exponent 1 at a vanishing base: the second derivative vanishes
  const H3 linear = pow( variable( 0.0, 1 ), 1 );
  checkClose( linear.dx( 1 ), 1.0, "d x^1/dx at 0" );
  checkClose( linear.dxx( 1, 1 ), 0.0, "d2 x^1/dx2 at 0" );

  // seeding a number again makes it the new variable alone
  H3 x = variable( 0.7, 0 ) * variable( 1.3, 1 );
  x.seed( 2 );
  checkClose( x.dx( 0 ) + x.dx( 1 ) + x.dxx( 0, 1 ) + x.dxx( 1, 0 ), 0.0, "reseeded: old derivatives" );
  checkClose( x.dx( 2 ), 1.0, "reseeded: new derivative" );
  throwExceptionOnFailure( x.mask() == 4u, Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " reseeded mask" );
}

void testSecondOrderAgainstHyperDual()
{
  // the full Hessian of a composite function of 3 variables against nested hyper-dual numbers, one pair at a time
  using Marmot::Math::Differentiation::Automatic::ScalarMode::HyperDual;
  const std::array< double, 3 > v = { 0.7, 1.3, 0.4 };
  const H3                      f = composite( variable( v[0], 0 ), variable( v[1], 1 ), variable( v[2], 2 ) );
  for ( std::size_t i = 0; i < 3; i++ )
    for ( std::size_t j = 0; j < 3; j++ ) {
      std::array< HyperDual< 2 >, 3 > x;
      for ( std::size_t k = 0; k < 3; k++ )
        x[k] = v[k];
      x[i].grad()            = 1.0;
      x[j].val().grad()      = 1.0;
      const HyperDual< 2 > g = composite( x[0], x[1], x[2] );
      checkClose( f.val(), g.val().val(), "composite value" );
      checkClose( f.dx( i ), g.grad().val(), "composite gradient" );
      checkClose( f.dxx( i, j ), g.grad().grad(), "composite Hessian" );
    }
}

namespace {

  using T3 = Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< 3 >;

  T3 thirdOrderVariable( double value, std::size_t direction )
  {
    T3 x = value;
    x.seed( direction );
    return x;
  }

} // namespace

void testThirdOrderLayout()
{
  // the sorted triples are stored in the order i, then j, then k, and every permutation of a triple is the same entry
  using T5      = Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< 5 >;
  std::size_t q = 0;
  for ( std::size_t i = 0; i < 5; i++ )
    for ( std::size_t j = i; j < 5; j++ ) {
      throwExceptionOnFailure( T5::thirdRowStart( i, j ) == q,
                               Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " row start" );
      for ( std::size_t k = j; k < 5; k++, q++ )
        throwExceptionOnFailure( T5::thirdIndex( i, j, k ) == q && T5::thirdIndex( i, k, j ) == q &&
                                   T5::thirdIndex( j, i, k ) == q && T5::thirdIndex( j, k, i ) == q &&
                                   T5::thirdIndex( k, i, j ) == q && T5::thirdIndex( k, j, i ) == q,
                                 Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " index of a triple" );
    }
  throwExceptionOnFailure( q == T5::nThird && q == 35,
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " number of triples" );
}

void testThirdOrderElementaryFunctions()
{
  // f( x ) with x the variable 1 of 3: value, f', f'' in ( 1, 1 ), f''' in ( 1, 1, 1 ), everything else zero
  const double x = 0.6;
  struct Case {
    std::string                      name;
    std::function< T3( const T3& ) > f;
    double                           f0, f1, f2, f3;
  };
  const double              t = std::tan( x ), th = std::tanh( x ), ln2 = std::log( 2 );
  const std::vector< Case > cases = {
    { "sqrt",
      []( const T3& a ) { return sqrt( a ); },
      std::sqrt( x ),
      0.5 / std::sqrt( x ),
      -0.25 * std::pow( x, -1.5 ),
      0.375 * std::pow( x, -2.5 ) },
    { "cbrt",
      []( const T3& a ) { return cbrt( a ); },
      std::cbrt( x ),
      std::pow( x, -2. / 3 ) / 3,
      -2. / 9 * std::pow( x, -5. / 3 ),
      10. / 27 * std::pow( x, -8. / 3 ) },
    { "exp", []( const T3& a ) { return exp( a ); }, std::exp( x ), std::exp( x ), std::exp( x ), std::exp( x ) },
    { "log", []( const T3& a ) { return log( a ); }, std::log( x ), 1 / x, -1 / ( x * x ), 2 / ( x * x * x ) },
    { "log10",
      []( const T3& a ) { return log10( a ); },
      std::log10( x ),
      1 / ( x * std::numbers::ln10 ),
      -1 / ( x * x * std::numbers::ln10 ),
      2 / ( x * x * x * std::numbers::ln10 ) },
    { "pow( x, 2.5 )",
      []( const T3& a ) { return pow( a, 2.5 ); },
      std::pow( x, 2.5 ),
      2.5 * std::pow( x, 1.5 ),
      3.75 * std::pow( x, 0.5 ),
      1.875 * std::pow( x, -0.5 ) },
    { "pow( x, 3 )", []( const T3& a ) { return pow( a, 3 ); }, x * x * x, 3 * x * x, 6 * x, 6 },
    { "pow( x, 2 )", []( const T3& a ) { return pow( a, 2 ); }, x * x, 2 * x, 2, 0 },
    { "pow( x, 0 )", []( const T3& a ) { return pow( a, 0 ); }, 1, 0, 0, 0 },
    { "pow( 2, x )",
      []( const T3& a ) { return pow( 2., a ); },
      std::pow( 2, x ),
      std::pow( 2, x ) * ln2,
      std::pow( 2, x ) * ln2 * ln2,
      std::pow( 2, x ) * ln2 * ln2 * ln2 },
    { "sin", []( const T3& a ) { return sin( a ); }, std::sin( x ), std::cos( x ), -std::sin( x ), -std::cos( x ) },
    { "cos", []( const T3& a ) { return cos( a ); }, std::cos( x ), -std::sin( x ), -std::cos( x ), std::sin( x ) },
    { "tan",
      []( const T3& a ) { return tan( a ); },
      t,
      1 + t * t,
      2 * t * ( 1 + t * t ),
      ( 2 + 4 * std::pow( std::sin( x ), 2 ) ) / std::pow( std::cos( x ), 4 ) },
    { "asin",
      []( const T3& a ) { return asin( a ); },
      std::asin( x ),
      std::pow( 1 - x * x, -0.5 ),
      x * std::pow( 1 - x * x, -1.5 ),
      ( 1 + 2 * x * x ) * std::pow( 1 - x * x, -2.5 ) },
    { "acos",
      []( const T3& a ) { return acos( a ); },
      std::acos( x ),
      -std::pow( 1 - x * x, -0.5 ),
      -x * std::pow( 1 - x * x, -1.5 ),
      -( 1 + 2 * x * x ) * std::pow( 1 - x * x, -2.5 ) },
    { "atan",
      []( const T3& a ) { return atan( a ); },
      std::atan( x ),
      1 / ( 1 + x * x ),
      -2 * x / std::pow( 1 + x * x, 2 ),
      ( 6 * x * x - 2 ) / std::pow( 1 + x * x, 3 ) },
    { "sinh", []( const T3& a ) { return sinh( a ); }, std::sinh( x ), std::cosh( x ), std::sinh( x ), std::cosh( x ) },
    { "cosh", []( const T3& a ) { return cosh( a ); }, std::cosh( x ), std::sinh( x ), std::cosh( x ), std::sinh( x ) },
    { "tanh",
      []( const T3& a ) { return tanh( a ); },
      th,
      1 - th * th,
      -2 * th * ( 1 - th * th ),
      ( 4 * std::pow( std::sinh( x ), 2 ) - 2 ) / std::pow( std::cosh( x ), 4 ) },
    { "abs( -x )", []( const T3& a ) { return abs( -a ); }, x, 1, 0, 0 },
    { "1 / x", []( const T3& a ) { return 1. / a; }, 1 / x, -1 / ( x * x ), 2 / ( x * x * x ), -6 / ( x * x * x * x ) },
    { "2 - 3 x", []( const T3& a ) { return 2. - a * 3.; }, 2 - 3 * x, -3, 0, 0 },
  };
  for ( const auto& c : cases ) {
    const T3 y = c.f( thirdOrderVariable( x, 1 ) );
    checkClose( y.val(), c.f0, c.name + " value" );
    checkClose( y.dx( 1 ), c.f1, c.name + " first derivative" );
    checkClose( y.dxx( 1, 1 ), c.f2, c.name + " second derivative" );
    checkClose( y.dxxx( 1, 1, 1 ), c.f3, c.name + " third derivative" );
    for ( std::size_t i = 0; i < 3; i++ )
      for ( std::size_t j = 0; j < 3; j++ ) {
        if ( i != 1 || j != 1 )
          checkClose( y.dxx( i, j ), 0.0, c.name + " other second derivatives" );
        for ( std::size_t k = 0; k < 3; k++ )
          if ( i != 1 || j != 1 || k != 1 )
            checkClose( y.dxxx( i, j, k ), 0.0, c.name + " other third derivatives" );
      }
    throwExceptionOnFailure( y.dx( 0 ) == 0.0 && y.dx( 2 ) == 0.0,
                             Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << c.name << " other derivatives" );
  }
}

void testThirdOrderBinaryOperations()
{
  // x and y are the variables 0 and 2 of 3, so the operands depend on different directions
  const double x = 0.7, y = 1.3;
  const T3     x_ = thirdOrderVariable( x, 0 ), y_ = thirdOrderVariable( y, 2 );

  const T3 p = x_ * y_;
  checkClose( p.val(), x * y, "x y" );
  checkClose( p.dx( 0 ), y, "d( x y )/dx" );
  checkClose( p.dx( 2 ), x, "d( x y )/dy" );
  checkClose( p.dxx( 0, 2 ), 1, "d2( x y )/dxdy" );
  checkClose( p.dxx( 2, 0 ), 1, "d2( x y )/dydx" );
  for ( std::size_t i = 0; i < 3; i++ )
    for ( std::size_t j = 0; j < 3; j++ )
      for ( std::size_t k = 0; k < 3; k++ )
        checkClose( p.dxxx( i, j, k ), 0, "d3( x y )" );

  const T3 xxy = x_ * x_ * y_;
  checkClose( xxy.dxxx( 0, 0, 2 ), 2, "d3( x x y )/dx2dy" );
  checkClose( xxy.dxxx( 0, 0, 0 ) + xxy.dxxx( 0, 2, 2 ) + xxy.dxxx( 2, 2, 2 ), 0, "d3( x x y ), other" );

  const T3 q = x_ / y_;
  checkClose( q.dx( 2 ), -x / ( y * y ), "d( x / y )/dy" );
  checkClose( q.dxx( 0, 2 ), -1 / ( y * y ), "d2( x / y )/dxdy" );
  checkClose( q.dxx( 2, 2 ), 2 * x / ( y * y * y ), "d2( x / y )/dy2" );
  checkClose( q.dxxx( 0, 2, 2 ), 2 / ( y * y * y ), "d3( x / y )/dxdy2" );
  checkClose( q.dxxx( 2, 0, 2 ), 2 / ( y * y * y ), "d3( x / y )/dydxdy" );
  checkClose( q.dxxx( 2, 2, 0 ), 2 / ( y * y * y ), "d3( x / y )/dy2dx" );
  checkClose( q.dxxx( 2, 2, 2 ), -6 * x / ( y * y * y * y ), "d3( x / y )/dy3" );
  checkClose( q.dxxx( 0, 0, 2 ), 0, "d3( x / y )/dx2dy" );

  const T3 w = pow( x_, y_ );
  checkClose( w.dxx( 0, 2 ), std::pow( x, y - 1 ) * ( 1 + y * std::log( x ) ), "d2 x^y/dxdy" );
  checkClose( w.dxxx( 0, 0, 2 ), std::pow( x, y - 2 ) * ( 2 * y - 1 + y * ( y - 1 ) * std::log( x ) ), "d3 x^y/dx2dy" );
  checkClose( w.dxxx( 2, 2, 2 ), std::pow( x, y ) * std::pow( std::log( x ), 3 ), "d3 x^y/dy3" );

  const T3     a  = atan2( y_, x_ );
  const double r2 = x * x + y * y;
  checkClose( a.val(), std::atan2( y, x ), "atan2" );
  checkClose( a.dx( 0 ), -y / r2, "d atan2/dx" );
  checkClose( a.dxx( 0, 2 ), ( y * y - x * x ) / ( r2 * r2 ), "d2 atan2/dxdy" );
  checkClose( a.dxxx( 0, 2, 2 ), ( 6 * x * x * y - 2 * y * y * y ) / ( r2 * r2 * r2 ), "d3 atan2/dxdy2" );

  T3 s = x_;
  s += y_;
  s -= 2. * x_;
  s *= y_;
  s /= x_;
  // ( y - x ) y / x
  checkClose( s.val(), ( y - x ) * y / x, "compound assignment" );
  checkClose( s.dxx( 0, 2 ), -2 * y / ( x * x ), "compound assignment d2/dxdy" );
  checkClose( s.dxxx( 0, 2, 2 ), -2 / ( x * x ), "compound assignment d3/dxdy2" );
  checkClose( s.dxxx( 0, 0, 2 ), 4 * y / ( x * x * x ), "compound assignment d3/dx2dy" );
  checkClose( s.dxxx( 0, 0, 0 ), -6 * y * y / ( x * x * x * x ), "compound assignment d3/dx3" );
  s *= 2.;
  s /= 4.;
  checkClose( s.dxxx( 0, 0, 2 ), 2 * y / ( x * x * x ), "compound assignment with constants d3/dx2dy" );
  throwExceptionOnFailure( x_ < y_ && x_ == x && 0 < x_ && s != x_,
                           Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " comparisons" );
}

void testThirdOrderSpecialCases()
{
  // a vanishing base: pow without dividing by the base
  const T3 p = pow( thirdOrderVariable( 0.0, 1 ), 3.5 );
  checkClose( p.val(), 0.0, "0^3.5" );
  checkClose( p.dx( 1 ), 0.0, "d x^3.5/dx at 0" );
  checkClose( p.dxx( 1, 1 ), 0.0, "d2 x^3.5/dx2 at 0" );
  checkClose( p.dxxx( 1, 1, 1 ), 0.0, "d3 x^3.5/dx3 at 0" );

  // integer exponents at a vanishing base: the derivatives of an order above the exponent vanish
  const T3 linear = pow( thirdOrderVariable( 0.0, 1 ), 1 );
  checkClose( linear.dx( 1 ), 1.0, "d x^1/dx at 0" );
  checkClose( linear.dxx( 1, 1 ), 0.0, "d2 x^1/dx2 at 0" );
  checkClose( linear.dxxx( 1, 1, 1 ), 0.0, "d3 x^1/dx3 at 0" );
  const T3 cubic = pow( thirdOrderVariable( 0.0, 1 ), 3 );
  checkClose( cubic.dx( 1 ) + cubic.dxx( 1, 1 ), 0.0, "d x^3/dx, d2 x^3/dx2 at 0" );
  checkClose( cubic.dxxx( 1, 1, 1 ), 6.0, "d3 x^3/dx3 at 0" );

  // seeding a number again makes it the new variable alone
  T3 x = thirdOrderVariable( 0.7, 0 ) * thirdOrderVariable( 1.3, 1 ) * thirdOrderVariable( 1.3, 1 );
  checkClose( x.dxxx( 0, 1, 1 ), 2.0, "before reseeding" );
  x.seed( 2 );
  checkClose( x.dx( 0 ) + x.dx( 1 ) + x.dxx( 0, 1 ) + x.dxx( 1, 0 ) + x.dxxx( 0, 1, 1 ),
              0.0,
              "reseeded: old derivatives" );
  checkClose( x.dx( 2 ), 1.0, "reseeded: new derivative" );
  throwExceptionOnFailure( x.mask() == 4u, Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " reseeded mask" );
}

void testThirdOrderAgainstHyperDual()
{
  // all third derivatives of a composite function of 3 variables against nested hyper-dual numbers, one triple at a
  // time
  using Marmot::Math::Differentiation::Automatic::ScalarMode::HyperDual;
  const std::array< double, 3 > v = { 0.7, 1.3, 0.4 };
  const T3 f = composite( thirdOrderVariable( v[0], 0 ), thirdOrderVariable( v[1], 1 ), thirdOrderVariable( v[2], 2 ) );
  for ( std::size_t i = 0; i < 3; i++ )
    for ( std::size_t j = 0; j < 3; j++ )
      for ( std::size_t k = 0; k < 3; k++ ) {
        std::array< HyperDual< 3 >, 3 > x;
        for ( std::size_t l = 0; l < 3; l++ )
          x[l] = v[l];
        x[i].grad()             = 1.0;
        x[j].val().grad()       = 1.0;
        x[k].val().val().grad() = 1.0;
        const HyperDual< 3 > g  = composite( x[0], x[1], x[2] );
        checkClose( f.val(), g.val().val().val(), "composite value" );
        checkClose( f.dx( i ), g.grad().val().val(), "composite gradient" );
        checkClose( f.dxx( i, j ), g.grad().grad().val(), "composite Hessian" );
        checkClose( f.dxxx( i, j, k ), g.grad().grad().grad(), "composite third derivatives" );
      }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testArithmetic,
                                                       testAssignmentAndAliasing,
                                                       testElementaryFunctions,
                                                       testEigen,
                                                       testFastor,
                                                       testSecondOrderElementaryFunctions,
                                                       testSecondOrderBinaryOperations,
                                                       testSecondOrderSpecialCases,
                                                       testSecondOrderAgainstHyperDual,
                                                       testThirdOrderLayout,
                                                       testThirdOrderElementaryFunctions,
                                                       testThirdOrderBinaryOperations,
                                                       testThirdOrderSpecialCases,
                                                       testThirdOrderAgainstHyperDual };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
