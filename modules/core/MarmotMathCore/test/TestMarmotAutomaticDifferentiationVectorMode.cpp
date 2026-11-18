#include "Marmot/MarmotAutomaticDifferentiationVectorMode.h"
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

int main()
{
  auto tests = std::vector< std::function< void() > >{ testArithmetic,
                                                       testAssignmentAndAliasing,
                                                       testElementaryFunctions,
                                                       testEigen,
                                                       testFastor };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
