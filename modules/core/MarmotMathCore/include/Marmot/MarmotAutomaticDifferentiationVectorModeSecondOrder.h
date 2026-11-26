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

/** @file MarmotAutomaticDifferentiationVectorModeSecondOrder.h
 *  @brief Marmot's own second-order vector-mode forward automatic differentiation: value, gradient and Hessian with
 *  respect to nDirections variables in one evaluation.
 *
 *  VectorMode::Dual2nd< nDirections > carries a value, its gradient and the upper triangle of its (symmetric) Hessian.
 *  It is tailored to SecondOrder::d2f_dT2, which then evaluates the differentiated function once instead of once per
 *  pair of directions. Every operation applies the second-order chain rule: for y = f( a ),
 *  grad y = f' grad a and hess y = f' hess a + f'' grad a grad a^T, and for products the product rule. Elementary
 *  functions are evaluated once per operation, at the value.
 *
 *  Each number tracks which directions it depends on (a bit mask). Operations only update the Hessian rows of these
 *  directions, so the many operations on the seeded inputs themselves (e.g. the products of a determinant) cost a
 *  fraction of a full Hessian update.
 *
 *  It is a plain (eagerly evaluated) number type, so it may be stored in `auto`. Use it through the backend-neutral
 *  vocabulary of MarmotAutomaticDifferentiationTypes.h (HessianDual, seedDirection, partialDerivative,
 *  secondPartialDerivative).
 */

#pragma once
#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include <Eigen/Core>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <ostream>
#include <type_traits>
#include <utility>

// The operations are marked [[gnu::always_inline]]: without it, clang leaves them as calls, which makes a Hessian
// about 1.45 times slower. MSVC ignores the attribute.
namespace Marmot::Math::Differentiation::Automatic::VectorMode {

  template < std::size_t nDirections_ >
  class Dual2nd;

  /** @brief True if T is a second-order vector-mode number
   *  @tparam T type to check
   */
  template < typename T >
  inline constexpr bool isDual2nd = false;

  /** @brief True for second-order vector-mode numbers
   *  @tparam nDirections_ number of directions
   */
  template < std::size_t nDirections_ >
  inline constexpr bool isDual2nd< Dual2nd< nDirections_ > > = true;

  /** @brief Second-order vector-mode number: value, gradient and Hessian with respect to nDirections_ variables
   *  @tparam nDirections_ number of independent variables (directions)
   *
   *  The Hessian is stored as its upper triangle, row by row. The gradient and the Hessian precede the value, so that
   *  pairs of them are 16-byte aligned for SIMD.
   *
   *  The mask holds the directions the number depends on; the gradient and Hessian entries of all other directions
   *  are zero.
   */
  template < std::size_t nDirections_ >
  class Dual2nd {
  public:
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = nDirections_;
    /** @brief Number of stored Hessian entries */
    static constexpr std::size_t nHessian = nDirections * ( nDirections + 1 ) / 2;
    static_assert( nDirections >= 1 && nDirections <= 64, "the direction mask has 64 bits" );
    /** @brief Mask of all directions */
    static constexpr std::uint64_t allDirections = nDirections == 64 ? ~std::uint64_t( 0 )
                                                                     : ( std::uint64_t( 1 ) << nDirections ) - 1;

    /** @brief Position of the Hessian entry ( i, j ) in the stored upper triangle
     *  @param i first direction
     *  @param j second direction
     *  @return index into hessian()
     */
    static constexpr std::size_t hessianIndex( std::size_t i, std::size_t j )
    {
      if ( i > j )
        std::swap( i, j );
      return i * nDirections - i * ( i - 1 ) / 2 + ( j - i );
    }

    /** @brief Calls f( i ) for every direction i in a mask
     *  @tparam F callable
     *  @param mask directions
     *  @param f function of the direction
     */
    template < typename F >
    [[gnu::always_inline]] static inline void forEachDirection( std::uint64_t mask, F&& f )
    {
      if ( mask == allDirections )
        for ( std::size_t i = 0; i < nDirections; i++ )
          f( i );
      else
        for ( ; mask; mask &= mask - 1 )
          f( static_cast< std::size_t >( std::countr_zero( mask ) ) );
    }

    /** @brief Index of the diagonal entry ( i, i ) in the stored upper triangle, where row i starts
     *  @param i direction
     *  @return index into hessian()
     */
    static constexpr std::size_t rowStart( std::size_t i ) { return i * nDirections - i * ( i - 1 ) / 2; }

    /** @brief Calls f( i, k ) for every k of row i of the stored upper triangle, i.e., the entries ( i, j ), j >= i,
     *  for all directions i in a mask; the rows of the other directions do not change in the operations
     *  @tparam F callable
     *  @param mask directions
     *  @param f function of the row direction and the Hessian index; the column is j = i + k - rowStart( i )
     */
    template < typename F >
    [[gnu::always_inline]] static inline void forEachRow( std::uint64_t mask, F&& f )
    {
      forEachDirection( mask, [&]( std::size_t i ) {
        const std::size_t k0 = rowStart( i );
        for ( std::size_t j = i; j < nDirections; j++ )
          f( i, j, k0 + j - i );
      } );
    }

    /** @brief Tag of the constructor that leaves the gradient and the Hessian uninitialized */
    struct Uninitialized {};

    /** @brief Creates a zero */
    Dual2nd() : g_{}, h_{}, val_( 0.0 ), mask_( 0 ) {}

    /** @brief Creates a number with uninitialized gradient and Hessian, for the implementation of the operations,
     *  which must write all entries
     */
    explicit Dual2nd( Uninitialized ) : val_( 0.0 ), mask_( allDirections ) {}

    /** @brief Creates the result of an operation that depends on the directions in a mask: uninitialized if all
     *  directions, since the operation writes all entries then, and zero otherwise, since it writes only the
     *  entries of the rows in the mask
     *  @param mask directions of the result
     *  @return the number, with the mask set
     */
    [[gnu::always_inline]] static inline Dual2nd resultFor( std::uint64_t mask )
    {
      if ( mask == allDirections )
        return Dual2nd( Uninitialized{} );
      Dual2nd r;
      r.mask_ = mask;
      return r;
    }

    /** @brief Creates a constant, with zero derivatives
     *  @tparam S arithmetic type
     *  @param value value
     */
    template < ScalarMode::Constant S >
    Dual2nd( S value ) : g_{}, h_{}, val_( static_cast< double >( value ) ), mask_( 0 )
    {
    }

    /** @brief Makes this number the independent variable `direction` (gradient = unit vector, zero Hessian)
     *  @param direction index of the variable
     */
    void seed( std::size_t direction )
    {
      if ( mask_ != 0 ) {
        g_.fill( 0.0 );
        h_.fill( 0.0 );
      }
      g_[direction] = 1.0;
      mask_         = std::uint64_t( 1 ) << direction;
    }

    /** @brief Directions this number depends on
     *  @return bit mask
     */
    std::uint64_t mask() const { return mask_; }

    /** @brief Sets the directions this number depends on, for the implementation of the operations
     *  @param mask bit mask
     */
    void setMask( std::uint64_t mask ) { mask_ = mask; }

    /** @brief Value
     *  @return value
     */
    double val() const { return val_; }

    /** @brief Value
     *  @return reference to the value
     */
    double& val() { return val_; }

    /** @brief First partial derivative
     *  @param i direction
     *  @return partial derivative with respect to variable i
     */
    double dx( std::size_t i ) const { return g_[i]; }

    /** @brief Second partial derivative
     *  @param i first direction
     *  @param j second direction
     *  @return second partial derivative with respect to variables i and j
     */
    double dxx( std::size_t i, std::size_t j ) const { return h_[hessianIndex( i, j )]; }

    /** @brief Gradient
     *  @return the first partial derivatives
     */
    const std::array< double, nDirections >& gradient() const { return g_; }

    /** @brief Gradient, for the implementation of the operations, which keep the mask consistent
     *  @return reference to the first partial derivatives
     */
    std::array< double, nDirections >& gradient() { return g_; }

    /** @brief Hessian
     *  @return the upper triangle of the second partial derivatives, row by row
     */
    const std::array< double, nHessian >& hessian() const { return h_; }

    /** @brief Hessian, for the implementation of the operations, which keep the mask consistent
     *  @return reference to the upper triangle of the second partial derivatives, row by row
     */
    std::array< double, nHessian >& hessian() { return h_; }

    /** @brief Primal value
     *  @return value
     */
    explicit operator double() const { return val_; }

    /** @brief Adds a number
     *  @param b summand
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual2nd& operator+=( const Dual2nd& b )
    {
      combine( b, 1.0 );
      return *this;
    }

    /** @brief Subtracts a number
     *  @param b subtrahend
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual2nd& operator-=( const Dual2nd& b )
    {
      combine( b, -1.0 );
      return *this;
    }

    /** @brief Multiplies by a number
     *  @param b factor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual2nd& operator*=( const Dual2nd& b );

    /** @brief Divides by a number
     *  @param b divisor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual2nd& operator/=( const Dual2nd& b );

    /** @brief Adds a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual2nd& operator+=( S c )
    {
      val_ += static_cast< double >( c );
      return *this;
    }

    /** @brief Subtracts a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual2nd& operator-=( S c )
    {
      val_ -= static_cast< double >( c );
      return *this;
    }

    /** @brief Multiplies by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual2nd& operator*=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ *= c_;
      forEachDirection( mask_, [&]( std::size_t i ) { g_[i] *= c_; } );
      forEachRow( mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] *= c_; } );
      return *this;
    }

    /** @brief Divides by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual2nd& operator/=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ /= c_;
      forEachDirection( mask_, [&]( std::size_t i ) { g_[i] /= c_; } );
      forEachRow( mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] /= c_; } );
      return *this;
    }

  private:
    /** @brief this += sign b
     *  @param b number
     *  @param sign +1 or -1
     */
    [[gnu::always_inline]] inline void combine( const Dual2nd& b, double sign )
    {
      val_ += sign * b.val_;
      forEachDirection( b.mask_, [&]( std::size_t i ) { g_[i] += sign * b.g_[i]; } );
      forEachRow( b.mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] += sign * b.h_[k]; } );
      mask_ |= b.mask_;
    }

    /** @brief Gradient */
    std::array< double, nDirections > g_;
    /** @brief Upper triangle of the Hessian, row by row */
    std::array< double, nHessian > h_;
    /** @brief Value */
    double val_;
    /** @brief Directions this number depends on */
    std::uint64_t mask_;
  };

  /** @brief Applies a function, given by its value and first and second derivative at a.val(), to a number
   *  @tparam N number of directions
   *  @param a argument
   *  @param f0 f( a )
   *  @param f1 f'( a )
   *  @param f2 f''( a )
   *  @return f( a ): gradient f1 grad a, Hessian f1 hess a + f2 grad a grad a^T
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > chain( const Dual2nd< N >& a, double f0, double f1, double f2 )
  {
    using D        = Dual2nd< N >;
    D r            = D::resultFor( a.mask() );
    r.val()        = f0;
    const auto& ag = a.gradient();
    const auto& ah = a.hessian();
    auto&       rg = r.gradient();
    auto&       rh = r.hessian();
    D::forEachDirection( a.mask(), [&]( std::size_t i ) { rg[i] = f1 * ag[i]; } );
    D::forEachRow( a.mask(),
                   [&]( std::size_t i, std::size_t j, std::size_t k ) { rh[k] = f1 * ah[k] + f2 * ag[i] * ag[j]; } );
    return r;
  }

  /** @brief Product: grad = a grad b + b grad a, hess = a hess b + b hess a + grad a grad b^T + grad b grad a^T
   *  @tparam N number of directions
   *  @param a first factor
   *  @param b second factor
   *  @return a * b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator*( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    using D                = Dual2nd< N >;
    const double        av = a.val(), bv = b.val();
    const std::uint64_t m = a.mask() | b.mask();
    D                   r = D::resultFor( m );
    r.val()               = av * bv;
    const auto& ag        = a.gradient();
    const auto& bg        = b.gradient();
    const auto& ah        = a.hessian();
    const auto& bh        = b.hessian();
    auto&       rg        = r.gradient();
    auto&       rh        = r.hessian();
    // rows outside m stay zero; within a row, the columns outside m compute to zero
    D::forEachDirection( m, [&]( std::size_t i ) { rg[i] = av * bg[i] + bv * ag[i]; } );
    D::forEachRow( m, [&]( std::size_t i, std::size_t j, std::size_t k ) {
      rh[k] = av * bh[k] + bv * ah[k] + ag[i] * bg[j] + bg[i] * ag[j];
    } );
    return r;
  }

  template < std::size_t nDirections_ >
  Dual2nd< nDirections_ >& Dual2nd< nDirections_ >::operator*=( const Dual2nd& b )
  {
    return *this = *this * b;
  }

  /** @brief Reciprocal
   *  @tparam N number of directions
   *  @param a argument
   *  @return 1 / a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > reciprocal( const Dual2nd< N >& a )
  {
    const double r = 1.0 / a.val();
    return chain( a, r, -r * r, 2.0 * r * r * r );
  }

  /** @brief Quotient
   *  @tparam N number of directions
   *  @param a dividend
   *  @param b divisor
   *  @return a / b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator/( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return a * reciprocal( b );
  }

  template < std::size_t nDirections_ >
  Dual2nd< nDirections_ >& Dual2nd< nDirections_ >::operator/=( const Dual2nd& b )
  {
    return *this = *this / b;
  }

  /** @brief Linear combination ca a + cb b, writing every entry of the result once
   *  @tparam N number of directions
   *  @param a first number
   *  @param ca factor of a
   *  @param b second number
   *  @param cb factor of b
   *  @return ca a + cb b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > linearCombination( const Dual2nd< N >& a,
                                                                double              ca,
                                                                const Dual2nd< N >& b,
                                                                double              cb )
  {
    using D               = Dual2nd< N >;
    const std::uint64_t m = a.mask() | b.mask();
    D                   r = D::resultFor( m );
    r.val()               = ca * a.val() + cb * b.val();
    const auto& ag        = a.gradient();
    const auto& bg        = b.gradient();
    const auto& ah        = a.hessian();
    const auto& bh        = b.hessian();
    auto&       rg        = r.gradient();
    auto&       rh        = r.hessian();
    D::forEachDirection( m, [&]( std::size_t i ) { rg[i] = ca * ag[i] + cb * bg[i]; } );
    D::forEachRow( m, [&]( std::size_t, std::size_t, std::size_t k ) { rh[k] = ca * ah[k] + cb * bh[k]; } );
    return r;
  }

  /** @brief Product with a constant, writing every entry of the result once
   *  @tparam N number of directions
   *  @param a number
   *  @param c constant
   *  @return c a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > scaled( const Dual2nd< N >& a, double c )
  {
    using D        = Dual2nd< N >;
    D r            = D::resultFor( a.mask() );
    r.val()        = c * a.val();
    const auto& ag = a.gradient();
    const auto& ah = a.hessian();
    auto&       rg = r.gradient();
    auto&       rh = r.hessian();
    D::forEachDirection( a.mask(), [&]( std::size_t i ) { rg[i] = c * ag[i]; } );
    D::forEachRow( a.mask(), [&]( std::size_t, std::size_t, std::size_t k ) { rh[k] = c * ah[k]; } );
    return r;
  }

  /** @brief Unary plus
   *  @tparam N number of directions
   *  @param a operand
   *  @return a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator+( const Dual2nd< N >& a )
  {
    return a;
  }

  /** @brief Negation
   *  @tparam N number of directions
   *  @param a operand
   *  @return -a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator-( const Dual2nd< N >& a )
  {
    return scaled( a, -1.0 );
  }

  /** @brief Sum
   *  @tparam N number of directions
   *  @param a first summand
   *  @param b second summand
   *  @return a + b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator+( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return linearCombination( a, 1.0, b, 1.0 );
  }

  /** @brief Difference
   *  @tparam N number of directions
   *  @param a minuend
   *  @param b subtrahend
   *  @return a - b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator-( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return linearCombination( a, 1.0, b, -1.0 );
  }

  /** @brief Sum with a constant
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a number
   *  @param c constant
   *  @return a + c
   */
  template < std::size_t N, ScalarMode::Constant S >
  [[gnu::always_inline]] inline Dual2nd< N > operator+( Dual2nd< N > a, S c )
  {
    return a += c;
  }

  /** @brief Sum with a constant
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param c constant
   *  @param a number
   *  @return c + a
   */
  template < ScalarMode::Constant S, std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator+( S c, Dual2nd< N > a )
  {
    return a += c;
  }

  /** @brief Difference with a constant
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a number
   *  @param c constant
   *  @return a - c
   */
  template < std::size_t N, ScalarMode::Constant S >
  [[gnu::always_inline]] inline Dual2nd< N > operator-( Dual2nd< N > a, S c )
  {
    return a -= c;
  }

  /** @brief Difference of a constant and a number
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param c constant
   *  @param a number
   *  @return c - a
   */
  template < ScalarMode::Constant S, std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator-( S c, Dual2nd< N > a )
  {
    a *= -1.0;
    return a += c;
  }

  /** @brief Product with a constant
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a number
   *  @param c constant
   *  @return a * c
   */
  template < std::size_t N, ScalarMode::Constant S >
  [[gnu::always_inline]] inline Dual2nd< N > operator*( const Dual2nd< N >& a, S c )
  {
    return scaled( a, static_cast< double >( c ) );
  }

  /** @brief Product with a constant
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param c constant
   *  @param a number
   *  @return c * a
   */
  template < ScalarMode::Constant S, std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator*( S c, const Dual2nd< N >& a )
  {
    return scaled( a, static_cast< double >( c ) );
  }

  /** @brief Quotient by a constant
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a dividend
   *  @param c divisor
   *  @return a / c
   */
  template < std::size_t N, ScalarMode::Constant S >
  [[gnu::always_inline]] inline Dual2nd< N > operator/( Dual2nd< N > a, S c )
  {
    return a /= c;
  }

  /** @brief Quotient of a constant and a number
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param c dividend
   *  @param a divisor
   *  @return c / a
   */
  template < ScalarMode::Constant S, std::size_t N >
  [[gnu::always_inline]] inline Dual2nd< N > operator/( S c, const Dual2nd< N >& a )
  {
    return reciprocal( a ) * static_cast< double >( c );
  }

  /** @brief An operand pair of a comparison involving at least one second-order vector-mode number
   *  @tparam A first operand type
   *  @tparam B second operand type
   */
  template < typename A, typename B >
  concept Dual2ndComparable = ( isDual2nd< A > && (isDual2nd< B > || ScalarMode::Constant< B >)) ||
                              ( ScalarMode::Constant< A > && isDual2nd< B > );

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a == value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator==( const A& a, const B& b )
  {
    return static_cast< double >( a ) == static_cast< double >( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a != value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator!=( const A& a, const B& b )
  {
    return static_cast< double >( a ) != static_cast< double >( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a < value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator<( const A& a, const B& b )
  {
    return static_cast< double >( a ) < static_cast< double >( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a <= value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator<=( const A& a, const B& b )
  {
    return static_cast< double >( a ) <= static_cast< double >( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a > value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator>( const A& a, const B& b )
  {
    return static_cast< double >( a ) > static_cast< double >( b );
  }

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a >= value of b
   */
  template < typename A, typename B >
    requires Dual2ndComparable< A, B >
  bool operator>=( const A& a, const B& b )
  {
    return static_cast< double >( a ) >= static_cast< double >( b );
  }

  /** @brief Square root
   *  @tparam N number of directions
   *  @param a argument
   *  @return sqrt( a )
   */
  template < std::size_t N >
  Dual2nd< N > sqrt( const Dual2nd< N >& a )
  {
    const double s = std::sqrt( a.val() ), d = 0.5 / s;
    return chain( a, s, d, -0.5 * d / a.val() );
  }

  /** @brief Cube root
   *  @tparam N number of directions
   *  @param a argument
   *  @return cbrt( a )
   */
  template < std::size_t N >
  Dual2nd< N > cbrt( const Dual2nd< N >& a )
  {
    const double c = std::cbrt( a.val() ), d = 1.0 / ( 3.0 * c * c );
    return chain( a, c, d, -2.0 / 3.0 * d / a.val() );
  }

  /** @brief Exponential function
   *  @tparam N number of directions
   *  @param a argument
   *  @return exp( a )
   */
  template < std::size_t N >
  Dual2nd< N > exp( const Dual2nd< N >& a )
  {
    const double e = std::exp( a.val() );
    return chain( a, e, e, e );
  }

  /** @brief Natural logarithm
   *  @tparam N number of directions
   *  @param a argument
   *  @return log( a )
   */
  template < std::size_t N >
  Dual2nd< N > log( const Dual2nd< N >& a )
  {
    const double r = 1.0 / a.val();
    return chain( a, std::log( a.val() ), r, -r * r );
  }

  /** @brief Decimal logarithm
   *  @tparam N number of directions
   *  @param a argument
   *  @return log10( a )
   */
  template < std::size_t N >
  Dual2nd< N > log10( const Dual2nd< N >& a )
  {
    const double r = 1.0 / ( a.val() * std::numbers::ln10 );
    return chain( a, std::log10( a.val() ), r, -r / a.val() );
  }

  /** @brief Power with a constant exponent
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a base
   *  @param n exponent
   *  @return pow( a, n )
   */
  template < std::size_t N, ScalarMode::Constant S >
  Dual2nd< N > pow( const Dual2nd< N >& a, S n )
  {
    const double n_ = static_cast< double >( n ), x = a.val();
    if ( n_ == 2.0 )
      return a * a;
    if ( n_ == 0.0 )
      return Dual2nd< N >( 1.0 );
    if ( x != 0.0 ) {
      const double p = std::pow( x, n_ ), d = n_ * p / x;
      return chain( a, p, d, ( n_ - 1.0 ) * d / x );
    }
    // a vanishing base: the derivatives without dividing by the base; the second one of a linear function vanishes
    const double c2 = n_ * ( n_ - 1.0 );
    return chain( a, std::pow( x, n_ ), n_ * std::pow( x, n_ - 1.0 ), c2 == 0.0 ? 0.0 : c2 * std::pow( x, n_ - 2.0 ) );
  }

  /** @brief Power of a constant base
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param b base
   *  @param a exponent
   *  @return pow( b, a )
   */
  template < ScalarMode::Constant S, std::size_t N >
  Dual2nd< N > pow( S b, const Dual2nd< N >& a )
  {
    const double b_ = static_cast< double >( b ), p = std::pow( b_, a.val() );
    const double l = b_ == 0.0 ? 0.0 : std::log( b_ );
    return chain( a, p, p * l, p * l * l );
  }

  /** @brief Power of two numbers
   *  @tparam N number of directions
   *  @param a base
   *  @param b exponent
   *  @return pow( a, b ) = exp( b log( a ) )
   */
  template < std::size_t N >
  Dual2nd< N > pow( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return exp( b * log( a ) );
  }

  /** @brief Sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return sin( a )
   */
  template < std::size_t N >
  Dual2nd< N > sin( const Dual2nd< N >& a )
  {
    const double s = std::sin( a.val() );
    return chain( a, s, std::cos( a.val() ), -s );
  }

  /** @brief Cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return cos( a )
   */
  template < std::size_t N >
  Dual2nd< N > cos( const Dual2nd< N >& a )
  {
    const double c = std::cos( a.val() );
    return chain( a, c, -std::sin( a.val() ), -c );
  }

  /** @brief Tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return tan( a )
   */
  template < std::size_t N >
  Dual2nd< N > tan( const Dual2nd< N >& a )
  {
    const double t = std::tan( a.val() ), d = 1.0 + t * t;
    return chain( a, t, d, 2.0 * t * d );
  }

  /** @brief Arc sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return asin( a )
   */
  template < std::size_t N >
  Dual2nd< N > asin( const Dual2nd< N >& a )
  {
    const double x = a.val(), r = 1.0 / ( 1.0 - x * x ), d = std::sqrt( r );
    return chain( a, std::asin( x ), d, x * d * r );
  }

  /** @brief Arc cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return acos( a )
   */
  template < std::size_t N >
  Dual2nd< N > acos( const Dual2nd< N >& a )
  {
    const double x = a.val(), r = 1.0 / ( 1.0 - x * x ), d = std::sqrt( r );
    return chain( a, std::acos( x ), -d, -x * d * r );
  }

  /** @brief Arc tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return atan( a )
   */
  template < std::size_t N >
  Dual2nd< N > atan( const Dual2nd< N >& a )
  {
    const double x = a.val(), d = 1.0 / ( 1.0 + x * x );
    return chain( a, std::atan( x ), d, -2.0 * x * d * d );
  }

  /** @brief Arc tangent of y / x, using the signs of both to determine the quadrant
   *  @tparam N number of directions
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < std::size_t N >
  Dual2nd< N > atan2( const Dual2nd< N >& y, const Dual2nd< N >& x )
  {
    // atan2 differs from atan( y / x ) and from -atan( x / y ) by piecewise constants only; the quotient with the
    // larger denominator is used
    Dual2nd< N > r = std::abs( x.val() ) >= std::abs( y.val() ) ? atan( y / x ) : -atan( x / y );
    r.val()        = std::atan2( y.val(), x.val() );
    return r;
  }

  /** @brief Hyperbolic sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return sinh( a )
   */
  template < std::size_t N >
  Dual2nd< N > sinh( const Dual2nd< N >& a )
  {
    const double s = std::sinh( a.val() );
    return chain( a, s, std::cosh( a.val() ), s );
  }

  /** @brief Hyperbolic cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return cosh( a )
   */
  template < std::size_t N >
  Dual2nd< N > cosh( const Dual2nd< N >& a )
  {
    const double c = std::cosh( a.val() );
    return chain( a, c, std::sinh( a.val() ), c );
  }

  /** @brief Hyperbolic tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return tanh( a )
   */
  template < std::size_t N >
  Dual2nd< N > tanh( const Dual2nd< N >& a )
  {
    const double t = std::tanh( a.val() ), d = 1.0 - t * t;
    return chain( a, t, d, -2.0 * t * d );
  }

  /** @brief Absolute value; its derivatives at 0 are those of a
   *  @tparam N number of directions
   *  @param a argument
   *  @return |a|
   */
  template < std::size_t N >
  Dual2nd< N > abs( const Dual2nd< N >& a )
  {
    return a.val() < 0.0 ? -a : a;
  }

  /** @brief Absolute value; its derivatives at 0 are those of a
   *  @tparam N number of directions
   *  @param a argument
   *  @return |a|
   */
  template < std::size_t N >
  Dual2nd< N > fabs( const Dual2nd< N >& a )
  {
    return abs( a );
  }

  /** @brief Larger of two numbers; for equal values, a
   *  @tparam N number of directions
   *  @param a first number
   *  @param b second number
   *  @return the larger number
   */
  template < std::size_t N >
  Dual2nd< N > max( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return a.val() < b.val() ? b : a;
  }

  /** @brief Smaller of two numbers; for equal values, a
   *  @tparam N number of directions
   *  @param a first number
   *  @param b second number
   *  @return the smaller number
   */
  template < std::size_t N >
  Dual2nd< N > min( const Dual2nd< N >& a, const Dual2nd< N >& b )
  {
    return b.val() < a.val() ? b : a;
  }

  /** @brief True if the value is not a number
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isnan of the value
   */
  template < std::size_t N >
  bool isnan( const Dual2nd< N >& a )
  {
    return std::isnan( a.val() );
  }

  /** @brief True if the value is infinite
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isinf of the value
   */
  template < std::size_t N >
  bool isinf( const Dual2nd< N >& a )
  {
    return std::isinf( a.val() );
  }

  /** @brief True if the value is finite
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isfinite of the value
   */
  template < std::size_t N >
  bool isfinite( const Dual2nd< N >& a )
  {
    return std::isfinite( a.val() );
  }

  /** @brief Writes the value
   *  @tparam N number of directions
   *  @param out stream
   *  @param a number
   *  @return out
   */
  template < std::size_t N >
  std::ostream& operator<<( std::ostream& out, const Dual2nd< N >& a )
  {
    return out << a.val();
  }

} // namespace Marmot::Math::Differentiation::Automatic::VectorMode

namespace Eigen {
  /** @brief Eigen support for second-order vector-mode numbers
   *  @tparam N number of directions
   */
  template < std::size_t N >
  struct NumTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N > > : NumTraits< double > {
    /** @brief Real type */
    using Real = Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N >;
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

  /** @brief Lets Eigen mix second-order vector-mode numbers with double
   *  @tparam N number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t N, typename BinaryOp >
  struct ScalarBinaryOpTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N >, double, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N >;
  };

  /** @brief Lets Eigen mix double with second-order vector-mode numbers
   *  @tparam N number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t N, typename BinaryOp >
  struct ScalarBinaryOpTraits< double, Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N >, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual2nd< N >;
  };
} // namespace Eigen
