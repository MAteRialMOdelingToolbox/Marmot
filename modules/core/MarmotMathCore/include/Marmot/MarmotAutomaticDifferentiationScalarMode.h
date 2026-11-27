/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

/** @file MarmotAutomaticDifferentiationScalarMode.h
 *  @brief Marmot's own scalar-mode forward automatic differentiation: hyper-dual and univariate Taylor numbers.
 *
 *  - ScalarMode::Dual< T > is a dual number whose value and derivative are of type T, either double or a Dual itself.
 *    Nesting gives hyper-dual numbers of any order (ScalarMode::HyperDual< order >), which carry one independent
 *    direction per nesting level and all mixed derivatives. Level 1 is the outermost one. A Dual is a dense array of
 *    2^order doubles: the value node, followed by the derivative node.
 *  - ScalarMode::Taylor< order > is a univariate truncated Taylor polynomial, which carries the derivatives 0..order
 *    of one variable.
 *
 *  Both are plain (eagerly evaluated) number types without expression templates, so they may be stored in `auto` and
 *  returned from functions with deduced return types. Elementary functions of Taylor numbers evaluate the
 *  transcendental function once, at the value, and compose its Taylor series with the argument.
 *
 *  Use this header through the backend-neutral vocabulary of MarmotAutomaticDifferentiationTypes.h (dual, dual2nd,
 *  dual3rd, HyperDual, Taylor, seed, derivative, value, taylorDerivatives).
 */

#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <ostream>
#include <type_traits>

// The elementary operations are marked [[gnu::always_inline]]: for third order hyper-dual numbers their nesting exceeds
// the inlining heuristics of GCC and clang, which then leave them as calls. MSVC ignores the attribute.
namespace Marmot::Math::Differentiation::Automatic::ScalarMode {

  /** @brief An arithmetic constant (e.g. double or int), which carries no derivatives
   *  @tparam S type to check
   */
  template < typename S >
  concept Constant = std::is_arithmetic_v< S >;

  /** @brief Calls f( std::integral_constant< std::size_t, i >{} ) for i = begin, ..., end - 1, unrolled at compile
   *  time, so that the triangular loops of the Taylor arithmetic become straight-line code
   *  @tparam begin first index
   *  @tparam end one past the last index
   *  @tparam F callable
   *  @param f function of the index
   */
  template < std::size_t begin, std::size_t end, typename F >
  [[gnu::always_inline]] inline void staticFor( F&& f )
  {
    if constexpr ( begin < end ) {
      f( std::integral_constant< std::size_t, begin >{} );
      staticFor< begin + 1, end >( f );
    }
  }

  template < typename T >
  class Dual;

  template < std::size_t order >
  class Taylor;

  /** @brief True if T is a (hyper-)dual number of any order
   *  @tparam T type to check
   */
  template < typename T >
  inline constexpr bool isDual = false;

  /** @brief True for (hyper-)dual numbers
   *  @tparam T type of the value and derivative
   */
  template < typename T >
  inline constexpr bool isDual< Dual< T > > = true;

  /** @brief Number of nesting levels of a (hyper-)dual number, 0 for double
   *  @tparam T number type
   */
  template < typename T >
  inline constexpr std::size_t orderOf = 0;

  /** @brief Number of nesting levels of a (hyper-)dual number
   *  @tparam T type of the value and derivative
   */
  template < typename T >
  inline constexpr std::size_t orderOf< Dual< T > > = 1 + orderOf< T >;

  /** @brief Implementation of HyperDual
   *  @tparam order number of nesting levels
   */
  template < std::size_t order >
  struct HyperDualOfOrder {
    /** @brief The hyper-dual number type */
    using type = Dual< typename HyperDualOfOrder< order - 1 >::type >;
  };

  /** @brief Implementation of HyperDual: order 0 is double */
  template <>
  struct HyperDualOfOrder< 0 > {
    /** @brief The number type */
    using type = double;
  };

  /** @brief Hyper-dual number with the given number of nesting levels
   *  @tparam order number of nesting levels, i.e., the highest derivative order it carries
   */
  template < std::size_t order >
  using HyperDual = typename HyperDualOfOrder< order >::type;

  /** @brief Primal value of a (hyper-)dual number, or a double itself
   *  @tparam T number type
   *  @param x number
   *  @return the innermost value
   */
  template < typename T >
  double primalValue( const T& x )
  {
    if constexpr ( isDual< T > )
      return primalValue( x.val() );
    else
      return static_cast< double >( x );
  }

  /** @brief Dual number with value and derivative of type T
   *  @tparam T double, or a Dual for hyper-dual numbers
   */
  template < typename T >
  class Dual {
  public:
    /** @brief Type of the value and the derivative */
    using Inner = T;

    /** @brief Creates a zero */
    Dual() : val_( 0.0 ), grad_( 0.0 ) {}

    /** @brief Creates a constant, with zero derivatives
     *  @tparam S arithmetic type
     *  @param value value
     */
    template < Constant S >
    Dual( S value ) : val_( static_cast< double >( value ) ), grad_( 0.0 )
    {
    }

    /** @brief Creates a dual number from its value and derivative
     *  @param value value node
     *  @param grad derivative node
     */
    Dual( const T& value, const T& grad ) : val_( value ), grad_( grad ) {}

    /** @brief Primal value
     *  @return the innermost value
     */
    explicit operator double() const { return primalValue( *this ); }

    /** @brief Value node
     *  @return value node
     */
    const T& val() const { return val_; }

    /** @brief Value node
     *  @return reference to the value node
     */
    T& val() { return val_; }

    /** @brief Derivative node
     *  @return derivative node
     */
    const T& grad() const { return grad_; }

    /** @brief Derivative node
     *  @return reference to the derivative node
     */
    T& grad() { return grad_; }

    /** @brief Adds a number
     *  @param b summand
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual& operator+=( const Dual& b )
    {
      val_ += b.val_;
      grad_ += b.grad_;
      return *this;
    }

    /** @brief Subtracts a number
     *  @param b subtrahend
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual& operator-=( const Dual& b )
    {
      val_ -= b.val_;
      grad_ -= b.grad_;
      return *this;
    }

    /** @brief Multiplies by a number
     *  @param b factor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual& operator*=( const Dual& b )
    {
      // in place on every nesting level; b.grad_ is copied first, since b may be this number
      const T bGrad = b.grad_;
      grad_ *= b.val_;
      grad_ += val_ * bGrad;
      val_ *= b.val_;
      return *this;
    }

    /** @brief Divides by a number
     *  @param b divisor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual& operator/=( const Dual& b ) { return *this = *this / b; }

    /** @brief Adds a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    [[gnu::always_inline]] inline Dual& operator+=( S c )
    {
      val_ += static_cast< double >( c );
      return *this;
    }

    /** @brief Subtracts a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    [[gnu::always_inline]] inline Dual& operator-=( S c )
    {
      val_ -= static_cast< double >( c );
      return *this;
    }

    /** @brief Multiplies by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    [[gnu::always_inline]] inline Dual& operator*=( S c )
    {
      val_ *= static_cast< double >( c );
      grad_ *= static_cast< double >( c );
      return *this;
    }

    /** @brief Divides by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    [[gnu::always_inline]] inline Dual& operator/=( S c )
    {
      val_ /= static_cast< double >( c );
      grad_ /= static_cast< double >( c );
      return *this;
    }

  private:
    /** @brief Value node */
    T val_;
    /** @brief Derivative node */
    T grad_;
  };

  // ---------------------------------------------------------------------------------------------------------------
  // arithmetic of dual numbers
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Unary plus
   *  @tparam T inner type
   *  @param a operand
   *  @return a
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator+( const Dual< T >& a )
  {
    return a;
  }

  /** @brief Negation
   *  @tparam T inner type
   *  @param a operand
   *  @return -a
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator-( const Dual< T >& a )
  {
    return Dual< T >( -a.val(), -a.grad() );
  }

  /** @brief Sum
   *  @tparam T inner type
   *  @param a first summand
   *  @param b second summand
   *  @return a + b
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator+( const Dual< T >& a, const Dual< T >& b )
  {
    return Dual< T >( a.val() + b.val(), a.grad() + b.grad() );
  }

  /** @brief Sum with a constant
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dual number
   *  @param c constant
   *  @return a + c
   */
  template < typename T, Constant S >
  [[gnu::always_inline]] inline Dual< T > operator+( const Dual< T >& a, S c )
  {
    return Dual< T >( a.val() + static_cast< double >( c ), a.grad() );
  }

  /** @brief Sum with a constant
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param c constant
   *  @param a dual number
   *  @return c + a
   */
  template < Constant S, typename T >
  [[gnu::always_inline]] inline Dual< T > operator+( S c, const Dual< T >& a )
  {
    return Dual< T >( static_cast< double >( c ) + a.val(), a.grad() );
  }

  /** @brief Difference
   *  @tparam T inner type
   *  @param a minuend
   *  @param b subtrahend
   *  @return a - b
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator-( const Dual< T >& a, const Dual< T >& b )
  {
    return Dual< T >( a.val() - b.val(), a.grad() - b.grad() );
  }

  /** @brief Difference with a constant
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dual number
   *  @param c constant
   *  @return a - c
   */
  template < typename T, Constant S >
  [[gnu::always_inline]] inline Dual< T > operator-( const Dual< T >& a, S c )
  {
    return Dual< T >( a.val() - static_cast< double >( c ), a.grad() );
  }

  /** @brief Difference of a constant and a dual number
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param c constant
   *  @param a dual number
   *  @return c - a
   */
  template < Constant S, typename T >
  [[gnu::always_inline]] inline Dual< T > operator-( S c, const Dual< T >& a )
  {
    return Dual< T >( static_cast< double >( c ) - a.val(), -a.grad() );
  }

  /** @brief Product
   *  @tparam T inner type
   *  @param a first factor
   *  @param b second factor
   *  @return a * b
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator*( const Dual< T >& a, const Dual< T >& b )
  {
    Dual< T > r = a;
    return r *= b;
  }

  /** @brief Product with a constant
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dual number
   *  @param c constant
   *  @return a * c
   */
  template < typename T, Constant S >
  [[gnu::always_inline]] inline Dual< T > operator*( const Dual< T >& a, S c )
  {
    const double c_ = static_cast< double >( c );
    return Dual< T >( a.val() * c_, a.grad() * c_ );
  }

  /** @brief Product with a constant
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param c constant
   *  @param a dual number
   *  @return c * a
   */
  template < Constant S, typename T >
  [[gnu::always_inline]] inline Dual< T > operator*( S c, const Dual< T >& a )
  {
    const double c_ = static_cast< double >( c );
    return Dual< T >( c_ * a.val(), c_ * a.grad() );
  }

  /** @brief Reciprocal of a (hyper-)dual number or a double, with a single floating-point division
   *  @tparam T number type
   *  @param x number
   *  @return 1 / x
   */
  template < typename T >
  [[gnu::always_inline]] inline T reciprocal( const T& x )
  {
    if constexpr ( isDual< T > ) {
      const auto r = reciprocal( x.val() );
      return T( r, -( r * r ) * x.grad() );
    }
    else
      return 1.0 / x;
  }

  /** @brief Quotient
   *  @tparam T inner type
   *  @param a dividend
   *  @param b divisor
   *  @return a / b
   */
  template < typename T >
  [[gnu::always_inline]] inline Dual< T > operator/( const Dual< T >& a, const Dual< T >& b )
  {
    if constexpr ( isDual< T > ) {
      const T inverse = reciprocal( b.val() );
      const T q       = a.val() * inverse;
      return Dual< T >( q, ( a.grad() - q * b.grad() ) * inverse );
    }
    else {
      const T q = a.val() / b.val();
      return Dual< T >( q, ( a.grad() - q * b.grad() ) / b.val() );
    }
  }

  /** @brief Quotient by a constant
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dividend
   *  @param c divisor
   *  @return a / c
   */
  template < typename T, Constant S >
  [[gnu::always_inline]] inline Dual< T > operator/( const Dual< T >& a, S c )
  {
    const double c_ = static_cast< double >( c );
    return Dual< T >( a.val() / c_, a.grad() / c_ );
  }

  /** @brief Quotient of a constant and a dual number
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param c dividend
   *  @param a divisor
   *  @return c / a
   */
  template < Constant S, typename T >
  [[gnu::always_inline]] inline Dual< T > operator/( S c, const Dual< T >& a )
  {
    return static_cast< double >( c ) * reciprocal( a );
  }

  // ---------------------------------------------------------------------------------------------------------------
  // comparisons (of the primal values)
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief An operand pair of a comparison involving at least one dual number
   *  @tparam A first operand type
   *  @tparam B second operand type
   */
  template < typename A, typename B >
  concept DualComparable = ( isDual< A > && (isDual< B > || Constant< B >)) || ( Constant< A > && isDual< B > );

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) == primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator==( const A& a, const B& b )
  {
    return primalValue( a ) == primalValue( b );
  }

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) != primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator!=( const A& a, const B& b )
  {
    return primalValue( a ) != primalValue( b );
  }

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) < primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator<( const A& a, const B& b )
  {
    return primalValue( a ) < primalValue( b );
  }

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) <= primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator<=( const A& a, const B& b )
  {
    return primalValue( a ) <= primalValue( b );
  }

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) > primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator>( const A& a, const B& b )
  {
    return primalValue( a ) > primalValue( b );
  }

  /** @brief Compares the primal values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return primalValue( a ) >= primalValue( b )
   */
  template < typename A, typename B >
    requires DualComparable< A, B >
  bool operator>=( const A& a, const B& b )
  {
    return primalValue( a ) >= primalValue( b );
  }

  /** @brief Absolute value; its derivative at 0 is taken as 1
   *  @tparam T inner type
   *  @param x argument
   *  @return |x|
   */
  template < typename T >
  Dual< T > abs( const Dual< T >& x )
  {
    return primalValue( x ) < 0.0 ? -x : x;
  }

  /** @brief Absolute value; its derivative at 0 is taken as 1
   *  @tparam T inner type
   *  @param x argument
   *  @return |x|
   */
  template < typename T >
  Dual< T > fabs( const Dual< T >& x )
  {
    return abs( x );
  }

  /** @brief Larger of two numbers; for equal values, a
   *  @tparam T inner type
   *  @param a first number
   *  @param b second number
   *  @return the larger number
   */
  template < typename T >
  Dual< T > max( const Dual< T >& a, const Dual< T >& b )
  {
    return primalValue( a ) < primalValue( b ) ? b : a;
  }

  /** @brief Larger of a dual number and a constant; for equal values, a
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dual number
   *  @param b constant
   *  @return the larger number
   */
  template < typename T, Constant S >
  Dual< T > max( const Dual< T >& a, S b )
  {
    return primalValue( a ) < static_cast< double >( b ) ? Dual< T >( b ) : a;
  }

  /** @brief Larger of a constant and a dual number; for equal values, b
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param a constant
   *  @param b dual number
   *  @return the larger number
   */
  template < Constant S, typename T >
  Dual< T > max( S a, const Dual< T >& b )
  {
    return max( b, a );
  }

  /** @brief Smaller of two numbers; for equal values, a
   *  @tparam T inner type
   *  @param a first number
   *  @param b second number
   *  @return the smaller number
   */
  template < typename T >
  Dual< T > min( const Dual< T >& a, const Dual< T >& b )
  {
    return primalValue( b ) < primalValue( a ) ? b : a;
  }

  /** @brief Smaller of a dual number and a constant; for equal values, a
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param a dual number
   *  @param b constant
   *  @return the smaller number
   */
  template < typename T, Constant S >
  Dual< T > min( const Dual< T >& a, S b )
  {
    return static_cast< double >( b ) < primalValue( a ) ? Dual< T >( b ) : a;
  }

  /** @brief Smaller of a constant and a dual number; for equal values, b
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param a constant
   *  @param b dual number
   *  @return the smaller number
   */
  template < Constant S, typename T >
  Dual< T > min( S a, const Dual< T >& b )
  {
    return min( b, a );
  }

  /** @brief True if the primal value is not a number
   *  @tparam T inner type
   *  @param x number
   *  @return std::isnan of the value
   */
  template < typename T >
  bool isnan( const Dual< T >& x )
  {
    return std::isnan( primalValue( x ) );
  }

  /** @brief True if the primal value is infinite
   *  @tparam T inner type
   *  @param x number
   *  @return std::isinf of the value
   */
  template < typename T >
  bool isinf( const Dual< T >& x )
  {
    return std::isinf( primalValue( x ) );
  }

  /** @brief True if the primal value is finite
   *  @tparam T inner type
   *  @param x number
   *  @return std::isfinite of the value
   */
  template < typename T >
  bool isfinite( const Dual< T >& x )
  {
    return std::isfinite( primalValue( x ) );
  }

  /** @brief Writes the primal value
   *  @tparam T inner type
   *  @param out stream
   *  @param x number
   *  @return out
   */
  template < typename T >
  std::ostream& operator<<( std::ostream& out, const Dual< T >& x )
  {
    return out << primalValue( x );
  }

  // ---------------------------------------------------------------------------------------------------------------
  // univariate Taylor numbers
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Univariate truncated Taylor polynomial x( t ) = sum_k c_k t^k, k = 0..order
   *  @tparam order highest derivative order
   *
   *  The coefficients are the normalized derivatives c_k = x^(k) / k!. A variable x0 is Taylor( x0 ) with c_1 = 1.
   */
  template < std::size_t order >
  class Taylor {
  public:
    /** @brief Creates a zero */
    Taylor() : c_{} {}

    /** @brief Creates a constant
     *  @tparam S arithmetic type
     *  @param value value
     */
    template < Constant S >
    Taylor( S value ) : c_{}
    {
      c_[0] = static_cast< double >( value );
    }

    /** @brief Value
     *  @return c_0
     */
    explicit operator double() const { return c_[0]; }

    /** @brief Normalized derivative
     *  @param k order
     *  @return c_k = x^(k) / k!
     */
    double operator[]( std::size_t k ) const { return c_[k]; }

    /** @brief Normalized derivative
     *  @param k order
     *  @return reference to c_k = x^(k) / k!
     */
    double& operator[]( std::size_t k ) { return c_[k]; }

    /** @brief Adds a Taylor number
     *  @param b summand
     *  @return this number
     */
    Taylor& operator+=( const Taylor& b )
    {
      for ( std::size_t k = 0; k <= order; k++ )
        c_[k] += b.c_[k];
      return *this;
    }

    /** @brief Subtracts a Taylor number
     *  @param b subtrahend
     *  @return this number
     */
    Taylor& operator-=( const Taylor& b )
    {
      for ( std::size_t k = 0; k <= order; k++ )
        c_[k] -= b.c_[k];
      return *this;
    }

    /** @brief Multiplies by a Taylor number
     *  @param b factor
     *  @return this number
     */
    Taylor& operator*=( const Taylor& b ) { return *this = *this * b; }

    /** @brief Divides by a Taylor number
     *  @param b divisor
     *  @return this number
     */
    Taylor& operator/=( const Taylor& b ) { return *this = *this / b; }

    /** @brief Adds a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Taylor& operator+=( S c )
    {
      c_[0] += static_cast< double >( c );
      return *this;
    }

    /** @brief Subtracts a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Taylor& operator-=( S c )
    {
      c_[0] -= static_cast< double >( c );
      return *this;
    }

    /** @brief Multiplies by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Taylor& operator*=( S c )
    {
      for ( double& ck : c_ )
        ck *= static_cast< double >( c );
      return *this;
    }

    /** @brief Divides by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < Constant S >
    Taylor& operator/=( S c )
    {
      for ( double& ck : c_ )
        ck /= static_cast< double >( c );
      return *this;
    }

  private:
    /** @brief Normalized derivatives */
    std::array< double, order + 1 > c_;
  };

  /** @brief True if T is a Taylor number
   *  @tparam T type to check
   */
  template < typename T >
  inline constexpr bool isTaylor = false;

  /** @brief True for Taylor numbers
   *  @tparam order highest derivative order
   */
  template < std::size_t order >
  inline constexpr bool isTaylor< Taylor< order > > = true;

  /** @brief Unary plus
   *  @tparam order highest derivative order
   *  @param a operand
   *  @return a
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator+( const Taylor< order >& a )
  {
    return a;
  }

  /** @brief Negation
   *  @tparam order highest derivative order
   *  @param a operand
   *  @return -a
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator-( const Taylor< order >& a )
  {
    Taylor< order > r = a;
    return r *= -1.0;
  }

  /** @brief Sum
   *  @tparam order highest derivative order
   *  @param a first summand
   *  @param b second summand
   *  @return a + b
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator+( Taylor< order > a, const Taylor< order >& b )
  {
    return a += b;
  }

  /** @brief Sum with a constant
   *  @tparam order highest derivative order
   *  @tparam S arithmetic type
   *  @param a Taylor number
   *  @param c constant
   *  @return a + c
   */
  template < std::size_t order, Constant S >
  [[gnu::always_inline]] inline Taylor< order > operator+( Taylor< order > a, S c )
  {
    return a += c;
  }

  /** @brief Sum with a constant
   *  @tparam S arithmetic type
   *  @tparam order highest derivative order
   *  @param c constant
   *  @param a Taylor number
   *  @return c + a
   */
  template < Constant S, std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator+( S c, Taylor< order > a )
  {
    return a += c;
  }

  /** @brief Difference
   *  @tparam order highest derivative order
   *  @param a minuend
   *  @param b subtrahend
   *  @return a - b
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator-( Taylor< order > a, const Taylor< order >& b )
  {
    return a -= b;
  }

  /** @brief Difference with a constant
   *  @tparam order highest derivative order
   *  @tparam S arithmetic type
   *  @param a Taylor number
   *  @param c constant
   *  @return a - c
   */
  template < std::size_t order, Constant S >
  [[gnu::always_inline]] inline Taylor< order > operator-( Taylor< order > a, S c )
  {
    return a -= c;
  }

  /** @brief Difference of a constant and a Taylor number
   *  @tparam S arithmetic type
   *  @tparam order highest derivative order
   *  @param c constant
   *  @param a Taylor number
   *  @return c - a
   */
  template < Constant S, std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator-( S c, const Taylor< order >& a )
  {
    return -a + c;
  }

  /** @brief Product (truncated Cauchy product)
   *  @tparam order highest derivative order
   *  @param a first factor
   *  @param b second factor
   *  @return a * b
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator*( const Taylor< order >& a, const Taylor< order >& b )
  {
    Taylor< order > r;
    staticFor< 0, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = 0.0;
      staticFor< 0, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s += a[j] * b[k - j];
      } );
      r[k] = s;
    } );
    return r;
  }

  /** @brief Product with a constant
   *  @tparam order highest derivative order
   *  @tparam S arithmetic type
   *  @param a Taylor number
   *  @param c constant
   *  @return a * c
   */
  template < std::size_t order, Constant S >
  [[gnu::always_inline]] inline Taylor< order > operator*( Taylor< order > a, S c )
  {
    return a *= c;
  }

  /** @brief Product with a constant
   *  @tparam S arithmetic type
   *  @tparam order highest derivative order
   *  @param c constant
   *  @param a Taylor number
   *  @return c * a
   */
  template < Constant S, std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator*( S c, Taylor< order > a )
  {
    return a *= c;
  }

  /** @brief Quotient
   *  @tparam order highest derivative order
   *  @param a dividend
   *  @param b divisor
   *  @return a / b
   */
  template < std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator/( const Taylor< order >& a, const Taylor< order >& b )
  {
    // a = q b, solved for q_k
    Taylor< order > q;
    const double    inverse = 1.0 / b[0];
    staticFor< 0, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = a[k];
      staticFor< 1, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s -= b[j] * q[k - j];
      } );
      q[k] = s * inverse;
    } );
    return q;
  }

  /** @brief Quotient by a constant
   *  @tparam order highest derivative order
   *  @tparam S arithmetic type
   *  @param a dividend
   *  @param c divisor
   *  @return a / c
   */
  template < std::size_t order, Constant S >
  [[gnu::always_inline]] inline Taylor< order > operator/( Taylor< order > a, S c )
  {
    return a /= c;
  }

  /** @brief Quotient of a constant and a Taylor number
   *  @tparam S arithmetic type
   *  @tparam order highest derivative order
   *  @param c dividend
   *  @param a divisor
   *  @return c / a
   */
  template < Constant S, std::size_t order >
  [[gnu::always_inline]] inline Taylor< order > operator/( S c, const Taylor< order >& a )
  {
    return Taylor< order >( c ) / a;
  }

  /** @brief An operand pair of a comparison involving at least one Taylor number
   *  @tparam A first operand type
   *  @tparam B second operand type
   */
  template < typename A, typename B >
  concept TaylorComparable = ( isTaylor< A > && (isTaylor< B > || Constant< B >)) || ( Constant< A > && isTaylor< B > );

  /** @brief Value of a Taylor number or a constant, for the comparison operators
   *  @tparam T Taylor or arithmetic type
   *  @param x number
   *  @return value of x
   */
  template < typename T >
  double taylorValue( const T& x )
  {
    return static_cast< double >( x );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a == value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator==( const A& a, const B& b )
  {
    return taylorValue( a ) == taylorValue( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a != value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator!=( const A& a, const B& b )
  {
    return taylorValue( a ) != taylorValue( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a < value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator<( const A& a, const B& b )
  {
    return taylorValue( a ) < taylorValue( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a <= value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator<=( const A& a, const B& b )
  {
    return taylorValue( a ) <= taylorValue( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a > value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator>( const A& a, const B& b )
  {
    return taylorValue( a ) > taylorValue( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a >= value of b
   */
  template < typename A, typename B >
    requires TaylorComparable< A, B >
  bool operator>=( const A& a, const B& b )
  {
    return taylorValue( a ) >= taylorValue( b );
  }

  /** @brief Composes a Taylor series f( x0 + d ) = sum_k f_k d^k with the argument x = x0 + d
   *  @tparam order highest derivative order
   *  @param x argument
   *  @param f normalized derivatives f_k = f^(k)( x0 ) / k! of the function at x0 = x[ 0 ]
   *  @return f( x )
   */
  template < std::size_t order >
  Taylor< order > compose( const Taylor< order >& x, const std::array< double, order + 1 >& f )
  {
    Taylor< order > d = x;
    d[0]              = 0.0;
    // Horner's scheme, from the highest order; d has no constant term, so d^k vanishes beyond k = order
    Taylor< order > r( f[order] );
    for ( std::size_t k = order; k-- > 0; ) {
      r = r * d;
      r[0] += f[k];
    }
    return r;
  }

  /** @brief Composes a function with known derivative with the argument, f( x ) = f( x0 ) + int f'( x ) dx
   *  @tparam order highest derivative order
   *  @param x argument
   *  @param f0 f( x[ 0 ] )
   *  @param df f'( x ), as a Taylor number
   *  @return f( x )
   */
  template < std::size_t order >
  Taylor< order > integrate( const Taylor< order >& x, double f0, const Taylor< order >& df )
  {
    // y' = f'( x ) x'  <=>  k y_k = sum_{j=1..k} j x_j df_{k-j}
    Taylor< order > y( f0 );
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = 0.0;
      staticFor< 1, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s += double( j ) * x[j] * df[k - j];
      } );
      y[k] = s * ( 1.0 / double( k ) );
    } );
    return y;
  }

  /** @brief Exponential function
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return exp( x )
   */
  template < std::size_t order >
  Taylor< order > exp( const Taylor< order >& x )
  {
    // y' = y x'
    Taylor< order > y( std::exp( x[0] ) );
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = 0.0;
      staticFor< 1, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s += double( j ) * x[j] * y[k - j];
      } );
      y[k] = s * ( 1.0 / double( k ) );
    } );
    return y;
  }

  /** @brief Natural logarithm
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return log( x )
   */
  template < std::size_t order >
  Taylor< order > log( const Taylor< order >& x )
  {
    // x y' = x'
    Taylor< order > y( std::log( x[0] ) );
    const double    inverse = 1.0 / x[0];
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = double( k ) * x[k];
      staticFor< 1, k >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s -= double( j ) * y[j] * x[k - j];
      } );
      y[k] = s * inverse * ( 1.0 / double( k ) );
    } );
    return y;
  }

  /** @brief Decimal logarithm
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return log10( x )
   */
  template < std::size_t order >
  Taylor< order > log10( const Taylor< order >& x )
  {
    Taylor< order > y = log( x ) / std::numbers::ln10;
    y[0]              = std::log10( x[0] );
    return y;
  }

  /** @brief Normalized derivatives of x^n at x0, with the value given separately
   *  @tparam order highest derivative order
   *  @param x0 point
   *  @param n exponent
   *  @param value x0^n
   *  @return binomial( n, k ) x0^(n - k), k = 0..order; for x0 = 0 exact for non-negative integer n
   */
  template < std::size_t order >
  std::array< double, order + 1 > powerSeries( double x0, double n, double value )
  {
    std::array< double, order + 1 > f{};
    f[0] = value;
    for ( std::size_t k = 1; k <= order; k++ )
      f[k] = x0 != 0.0 ? f[k - 1] * ( n - double( k - 1 ) ) / ( double( k ) * x0 ) : ( n == double( k ) ? 1.0 : 0.0 );
    return f;
  }

  /** @brief Power x^n with the value x0^n given separately (e.g. as sqrt( x0 ) for n = 1/2)
   *  @tparam order highest derivative order
   *  @param x base
   *  @param n exponent
   *  @param value x[ 0 ]^n
   *  @return x^n
   */
  template < std::size_t order >
  Taylor< order > power( const Taylor< order >& x, double n, double value )
  {
    if ( x[0] == 0.0 )
      return compose( x, powerSeries< order >( 0.0, n, value ) );
    // x y' = n y x'
    Taylor< order > y( value );
    const double    inverse = 1.0 / x[0];
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = 0.0;
      staticFor< 0, k >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s += ( n * double( k - j ) - double( j ) ) * x[k - j] * y[j];
      } );
      y[k] = s * inverse * ( 1.0 / double( k ) );
    } );
    return y;
  }

  /** @brief Power with a constant exponent
   *  @tparam order highest derivative order
   *  @tparam S arithmetic type
   *  @param x base
   *  @param n exponent
   *  @return pow( x, n )
   */
  template < std::size_t order, Constant S >
  Taylor< order > pow( const Taylor< order >& x, S n )
  {
    const double n_ = static_cast< double >( n );
    // small integer exponents by multiplication, which is faster and exact
    if ( n_ == 2.0 )
      return x * x;
    if ( n_ == 3.0 )
      return x * x * x;
    if ( n_ == 4.0 ) {
      const Taylor< order > x2 = x * x;
      return x2 * x2;
    }
    return power( x, n_, std::pow( x[0], n_ ) );
  }

  /** @brief Power of a constant base
   *  @tparam S arithmetic type
   *  @tparam order highest derivative order
   *  @param b base
   *  @param x exponent
   *  @return pow( b, x )
   */
  template < Constant S, std::size_t order >
  Taylor< order > pow( S b, const Taylor< order >& x )
  {
    const double    b_ = static_cast< double >( b );
    Taylor< order > y  = exp( x * std::log( b_ ) );
    y[0]               = std::pow( b_, x[0] );
    return y;
  }

  /** @brief Power of two Taylor numbers
   *  @tparam order highest derivative order
   *  @param x base
   *  @param y exponent
   *  @return pow( x, y )
   */
  template < std::size_t order >
  Taylor< order > pow( const Taylor< order >& x, const Taylor< order >& y )
  {
    return exp( y * log( x ) );
  }

  /** @brief Square root
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return sqrt( x )
   */
  template < std::size_t order >
  Taylor< order > sqrt( const Taylor< order >& x )
  {
    if ( x[0] == 0.0 )
      return compose( x, powerSeries< order >( 0.0, 0.5, 0.0 ) );
    // y^2 = x
    Taylor< order > y( std::sqrt( x[0] ) );
    const double    inverse = 0.5 / y[0];
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = x[k];
      staticFor< 1, k >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s -= y[j] * y[k - j];
      } );
      y[k] = s * inverse;
    } );
    return y;
  }

  /** @brief Cube root
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return cbrt( x )
   */
  template < std::size_t order >
  Taylor< order > cbrt( const Taylor< order >& x )
  {
    return power( x, 1.0 / 3.0, std::cbrt( x[0] ) );
  }

  /** @brief Joint recurrence of a pair of functions with s' = c x' and c' = sign s x' (sin/cos, sinh/cosh)
   *  @tparam order highest derivative order
   *  @param x argument
   *  @param s s( x[ 0 ] )
   *  @param c c( x[ 0 ] )
   *  @param sign -1 for sin/cos, +1 for sinh/cosh
   *  @return { s( x ), c( x ) }
   */
  template < std::size_t order >
  std::array< Taylor< order >, 2 > sineCosine( const Taylor< order >& x, double s, double c, double sign )
  {
    std::array< Taylor< order >, 2 > sc{ Taylor< order >( s ), Taylor< order >( c ) };
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k  = decltype( k_ )::value;
      double                ds = 0.0, dc = 0.0;
      staticFor< 1, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        ds += double( j ) * x[j] * sc[1][k - j];
        dc += double( j ) * x[j] * sc[0][k - j];
      } );
      sc[0][k] = ds * ( 1.0 / double( k ) );
      sc[1][k] = sign * dc * ( 1.0 / double( k ) );
    } );
    return sc;
  }

  /** @brief Sine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return sin( x )
   */
  template < std::size_t order >
  Taylor< order > sin( const Taylor< order >& x )
  {
    return sineCosine( x, std::sin( x[0] ), std::cos( x[0] ), -1.0 )[0];
  }

  /** @brief Cosine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return cos( x )
   */
  template < std::size_t order >
  Taylor< order > cos( const Taylor< order >& x )
  {
    return sineCosine( x, std::sin( x[0] ), std::cos( x[0] ), -1.0 )[1];
  }

  /** @brief Hyperbolic sine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return sinh( x )
   */
  template < std::size_t order >
  Taylor< order > sinh( const Taylor< order >& x )
  {
    return sineCosine( x, std::sinh( x[0] ), std::cosh( x[0] ), 1.0 )[0];
  }

  /** @brief Hyperbolic cosine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return cosh( x )
   */
  template < std::size_t order >
  Taylor< order > cosh( const Taylor< order >& x )
  {
    return sineCosine( x, std::sinh( x[0] ), std::cosh( x[0] ), 1.0 )[1];
  }

  /** @brief Function y with y' = ( 1 + sign y^2 ) x' (tan: sign = +1, tanh: sign = -1)
   *  @tparam order highest derivative order
   *  @param x argument
   *  @param y0 y( x[ 0 ] )
   *  @param sign +1 or -1
   *  @return y( x )
   */
  template < std::size_t order >
  Taylor< order > tangent( const Taylor< order >& x, double y0, double sign )
  {
    Taylor< order > y( y0 );
    Taylor< order > g; // 1 + sign y^2, completed up to order k - 1 before y_k is computed
    g[0] = 1.0 + sign * y0 * y0;
    staticFor< 1, order + 1 >( [&]( auto k_ ) {
      constexpr std::size_t k = decltype( k_ )::value;
      double                s = 0.0;
      staticFor< 1, k + 1 >( [&]( auto j_ ) {
        constexpr std::size_t j = decltype( j_ )::value;
        s += double( j ) * x[j] * g[k - j];
      } );
      y[k] = s * ( 1.0 / double( k ) );
      if constexpr ( k < order ) {
        double y2 = 0.0;
        staticFor< 0, k + 1 >( [&]( auto i_ ) {
          constexpr std::size_t i = decltype( i_ )::value;
          y2 += y[i] * y[k - i];
        } );
        g[k] = sign * y2;
      }
    } );
    return y;
  }

  /** @brief Tangent
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return tan( x )
   */
  template < std::size_t order >
  Taylor< order > tan( const Taylor< order >& x )
  {
    return tangent( x, std::tan( x[0] ), 1.0 );
  }

  /** @brief Hyperbolic tangent
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return tanh( x )
   */
  template < std::size_t order >
  Taylor< order > tanh( const Taylor< order >& x )
  {
    return tangent( x, std::tanh( x[0] ), -1.0 );
  }

  /** @brief Arc sine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return asin( x )
   */
  template < std::size_t order >
  Taylor< order > asin( const Taylor< order >& x )
  {
    return integrate( x, std::asin( x[0] ), 1.0 / sqrt( 1.0 - x * x ) );
  }

  /** @brief Arc cosine
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return acos( x )
   */
  template < std::size_t order >
  Taylor< order > acos( const Taylor< order >& x )
  {
    return integrate( x, std::acos( x[0] ), -1.0 / sqrt( 1.0 - x * x ) );
  }

  /** @brief Arc tangent
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return atan( x )
   */
  template < std::size_t order >
  Taylor< order > atan( const Taylor< order >& x )
  {
    return integrate( x, std::atan( x[0] ), 1.0 / ( 1.0 + x * x ) );
  }

  /** @brief Absolute value; its derivatives at 0 are those of x
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return |x|
   */
  template < std::size_t order >
  Taylor< order > abs( const Taylor< order >& x )
  {
    return x[0] < 0.0 ? -x : x;
  }

  /** @brief Absolute value; its derivatives at 0 are those of x
   *  @tparam order highest derivative order
   *  @param x argument
   *  @return |x|
   */
  template < std::size_t order >
  Taylor< order > fabs( const Taylor< order >& x )
  {
    return abs( x );
  }

  /** @brief Writes the value
   *  @tparam order highest derivative order
   *  @param out stream
   *  @param x number
   *  @return out
   */
  template < std::size_t order >
  std::ostream& operator<<( std::ostream& out, const Taylor< order >& x )
  {
    return out << x[0];
  }

  // ---------------------------------------------------------------------------------------------------------------
  // elementary functions of dual numbers: chain rule f( a + b eps ) = f( a ) + f'( a ) b eps, recursively, with
  // f'( a ) expressed by f( a ) where possible, so that each function calls its double counterpart once
  // ---------------------------------------------------------------------------------------------------------------

  /** @brief Square root
   *  @tparam T inner type
   *  @param x argument
   *  @return sqrt( x )
   */
  template < typename T >
  Dual< T > sqrt( const Dual< T >& x )
  {
    using std::sqrt;
    const T s = sqrt( x.val() );
    return Dual< T >( s, ( 0.5 * reciprocal( s ) ) * x.grad() );
  }

  /** @brief Cube root
   *  @tparam T inner type
   *  @param x argument
   *  @return cbrt( x )
   */
  template < typename T >
  Dual< T > cbrt( const Dual< T >& x )
  {
    using std::cbrt;
    const T c = cbrt( x.val() );
    return Dual< T >( c, ( ( 1.0 / 3.0 ) * reciprocal( c * c ) ) * x.grad() );
  }

  /** @brief Exponential function
   *  @tparam T inner type
   *  @param x argument
   *  @return exp( x )
   */
  template < typename T >
  Dual< T > exp( const Dual< T >& x )
  {
    using std::exp;
    const T e = exp( x.val() );
    return Dual< T >( e, e * x.grad() );
  }

  /** @brief Natural logarithm
   *  @tparam T inner type
   *  @param x argument
   *  @return log( x )
   */
  template < typename T >
  Dual< T > log( const Dual< T >& x )
  {
    using std::log;
    return Dual< T >( log( x.val() ), reciprocal( x.val() ) * x.grad() );
  }

  /** @brief Decimal logarithm
   *  @tparam T inner type
   *  @param x argument
   *  @return log10( x )
   */
  template < typename T >
  Dual< T > log10( const Dual< T >& x )
  {
    using std::log10;
    return Dual< T >( log10( x.val() ), ( ( 1.0 / std::numbers::ln10 ) * reciprocal( x.val() ) ) * x.grad() );
  }

  /** @brief Power with a constant exponent
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param x base
   *  @param n exponent
   *  @return pow( x, n )
   */
  template < typename T, Constant S >
  Dual< T > pow( const Dual< T >& x, S n )
  {
    using std::pow;
    const double n_ = static_cast< double >( n );
    // small integer exponents by multiplication, which is faster and exact
    if ( n_ == 2.0 )
      return x * x;
    if ( n_ == 3.0 )
      return x * x * x;
    if ( n_ == 4.0 ) {
      const Dual< T > x2 = x * x;
      return x2 * x2;
    }
    if ( n_ == 0.0 )
      return Dual< T >( 1.0 );
    const T p = pow( x.val(), n_ );
    if ( primalValue( x ) != 0.0 )
      return Dual< T >( p, ( n_ * p * reciprocal( x.val() ) ) * x.grad() );
    return Dual< T >( p, ( n_ * pow( x.val(), n_ - 1.0 ) ) * x.grad() );
  }

  /** @brief Power of a constant base
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param b base
   *  @param x exponent
   *  @return pow( b, x )
   */
  template < Constant S, typename T >
  Dual< T > pow( S b, const Dual< T >& x )
  {
    using std::pow;
    const double b_ = static_cast< double >( b );
    const T      p  = pow( b_, x.val() );
    if ( b_ == 0.0 )
      return Dual< T >( p, T( 0.0 ) );
    return Dual< T >( p, ( p * std::log( b_ ) ) * x.grad() );
  }

  /** @brief Power of two dual numbers
   *  @tparam T inner type
   *  @param x base
   *  @param y exponent
   *  @return pow( x, y )
   */
  template < typename T >
  Dual< T > pow( const Dual< T >& x, const Dual< T >& y )
  {
    using std::log;
    using std::pow;
    const T p = pow( x.val(), y.val() );
    if ( primalValue( x ) == 0.0 )
      return Dual< T >( p, y.val() * pow( x.val(), y.val() - 1.0 ) * x.grad() );
    return Dual< T >( p, p * ( y.val() * reciprocal( x.val() ) * x.grad() + log( x.val() ) * y.grad() ) );
  }

  /** @brief Sine and cosine of a (hyper-)dual number or a double, with one sin and one cos call
   *  @tparam T number type
   *  @param x argument
   *  @param sign -1 for sine and cosine, +1 for the hyperbolic sine and cosine
   *  @return { sin( x ), cos( x ) }, or { sinh( x ), cosh( x ) }
   */
  template < typename T >
  std::array< T, 2 > sineCosine( const T& x, double sign )
  {
    if constexpr ( isDual< T > ) {
      const auto sc = sineCosine( x.val(), sign );
      return { T( sc[0], sc[1] * x.grad() ), T( sc[1], sign * sc[0] * x.grad() ) };
    }
    else if ( sign < 0 )
      return { std::sin( x ), std::cos( x ) };
    else
      return { std::sinh( x ), std::cosh( x ) };
  }

  /** @brief Sine
   *  @tparam T inner type
   *  @param x argument
   *  @return sin( x )
   */
  template < typename T >
  Dual< T > sin( const Dual< T >& x )
  {
    return sineCosine( x, -1.0 )[0];
  }

  /** @brief Cosine
   *  @tparam T inner type
   *  @param x argument
   *  @return cos( x )
   */
  template < typename T >
  Dual< T > cos( const Dual< T >& x )
  {
    return sineCosine( x, -1.0 )[1];
  }

  /** @brief Tangent
   *  @tparam T inner type
   *  @param x argument
   *  @return tan( x )
   */
  template < typename T >
  Dual< T > tan( const Dual< T >& x )
  {
    using std::tan;
    const T t = tan( x.val() );
    return Dual< T >( t, ( 1.0 + t * t ) * x.grad() );
  }

  /** @brief Arc sine
   *  @tparam T inner type
   *  @param x argument
   *  @return asin( x )
   */
  template < typename T >
  Dual< T > asin( const Dual< T >& x )
  {
    using std::asin;
    using std::sqrt;
    return Dual< T >( asin( x.val() ), reciprocal( sqrt( 1.0 - x.val() * x.val() ) ) * x.grad() );
  }

  /** @brief Arc cosine
   *  @tparam T inner type
   *  @param x argument
   *  @return acos( x )
   */
  template < typename T >
  Dual< T > acos( const Dual< T >& x )
  {
    using std::acos;
    using std::sqrt;
    return Dual< T >( acos( x.val() ), -reciprocal( sqrt( 1.0 - x.val() * x.val() ) ) * x.grad() );
  }

  /** @brief Arc tangent
   *  @tparam T inner type
   *  @param x argument
   *  @return atan( x )
   */
  template < typename T >
  Dual< T > atan( const Dual< T >& x )
  {
    using std::atan;
    return Dual< T >( atan( x.val() ), reciprocal( 1.0 + x.val() * x.val() ) * x.grad() );
  }

  /** @brief Arc tangent of y / x, using the signs of both to determine the quadrant
   *  @tparam T inner type
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < typename T >
  Dual< T > atan2( const Dual< T >& y, const Dual< T >& x )
  {
    using std::atan2;
    const T r2 = x.val() * x.val() + y.val() * y.val();
    return Dual< T >( atan2( y.val(), x.val() ), ( x.val() * y.grad() - y.val() * x.grad() ) * reciprocal( r2 ) );
  }

  /** @brief Arc tangent of y / x with a constant x
   *  @tparam T inner type
   *  @tparam S arithmetic type
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < typename T, Constant S >
  Dual< T > atan2( const Dual< T >& y, S x )
  {
    return atan2( y, Dual< T >( x ) );
  }

  /** @brief Arc tangent of y / x with a constant y
   *  @tparam S arithmetic type
   *  @tparam T inner type
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < Constant S, typename T >
  Dual< T > atan2( S y, const Dual< T >& x )
  {
    return atan2( Dual< T >( y ), x );
  }

  /** @brief Hyperbolic sine
   *  @tparam T inner type
   *  @param x argument
   *  @return sinh( x )
   */
  template < typename T >
  Dual< T > sinh( const Dual< T >& x )
  {
    return sineCosine( x, 1.0 )[0];
  }

  /** @brief Hyperbolic cosine
   *  @tparam T inner type
   *  @param x argument
   *  @return cosh( x )
   */
  template < typename T >
  Dual< T > cosh( const Dual< T >& x )
  {
    return sineCosine( x, 1.0 )[1];
  }

  /** @brief Hyperbolic tangent
   *  @tparam T inner type
   *  @param x argument
   *  @return tanh( x )
   */
  template < typename T >
  Dual< T > tanh( const Dual< T >& x )
  {
    using std::tanh;
    const T t = tanh( x.val() );
    return Dual< T >( t, ( 1.0 - t * t ) * x.grad() );
  }

} // namespace Marmot::Math::Differentiation::Automatic::ScalarMode

namespace Eigen {
  /** @brief Eigen support for (hyper-)dual numbers
   *  @tparam T inner type
   */
  template < typename T >
  struct NumTraits< Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T > > : NumTraits< double > {
    /** @brief Real type */
    using Real = Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T >;
    /** @brief Non-integer type */
    using NonInteger = Real;
    /** @brief Literal type */
    using Literal = double;
    /** @brief Nested type */
    using Nested = Real;
    /** @brief Properties */
    enum {
      IsComplex             = 0,
      IsInteger             = 0,
      IsSigned              = 1,
      RequireInitialization = 1,
      ReadCost              = 1,
      AddCost               = 3,
      MulCost               = 3
    };
  };

  /** @brief Lets Eigen mix (hyper-)dual numbers with double
   *  @tparam T inner type
   *  @tparam BinaryOp Eigen binary functor
   */
  template < typename T, typename BinaryOp >
  struct ScalarBinaryOpTraits< Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T >, double, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T >;
  };

  /** @brief Lets Eigen mix double with (hyper-)dual numbers
   *  @tparam T inner type
   *  @tparam BinaryOp Eigen binary functor
   */
  template < typename T, typename BinaryOp >
  struct ScalarBinaryOpTraits< double, Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T >, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::ScalarMode::Dual< T >;
  };
} // namespace Eigen
