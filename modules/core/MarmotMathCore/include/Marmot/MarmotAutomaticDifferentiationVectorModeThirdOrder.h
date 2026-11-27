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

/** @file MarmotAutomaticDifferentiationVectorModeThirdOrder.h
 *  @brief Marmot's own third-order vector-mode forward automatic differentiation: value, gradient, Hessian and third
 *  derivatives with respect to nDirections variables in one evaluation.
 *
 *  VectorMode::Dual3rd< nDirections > carries a value, its gradient, the upper triangle of its (symmetric) Hessian and
 *  its third derivatives for the sorted triples of directions. It is tailored to ThirdOrder::d3f_dT3, which then
 *  evaluates the differentiated function once instead of once per triple of directions. Every operation applies the
 *  third-order chain rule: for y = f( a ),
 *  y_ijk = f' a_ijk + f'' ( a_i a_jk + a_j a_ik + a_k a_ij ) + f''' a_i a_j a_k, and for products the product rule.
 *  Elementary functions are evaluated once per operation, at the value.
 *
 *  Each number tracks which directions it depends on (a bit mask), like VectorMode::Dual2nd. Operations only update
 *  the Hessian rows and the rows of third derivatives of these directions.
 *
 *  It is a plain (eagerly evaluated) number type, so it may be stored in `auto`. Use it through the backend-neutral
 *  vocabulary of MarmotAutomaticDifferentiationTypes.h (ThirdDerivativeDual, seedDirection, partialDerivative,
 *  secondPartialDerivative, thirdPartialDerivative).
 */

#pragma once
#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorModeSecondOrder.h"
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

// The operations are marked [[gnu::always_inline]], like those of the second-order numbers. MSVC ignores the
// attribute.
namespace Marmot::Math::Differentiation::Automatic::VectorMode {

  template < std::size_t nDirections_ >
  class Dual3rd;

  /** @brief True if T is a third-order vector-mode number
   *  @tparam T type to check
   */
  template < typename T >
  inline constexpr bool isDual3rd = false;

  /** @brief True for third-order vector-mode numbers
   *  @tparam nDirections_ number of directions
   */
  template < std::size_t nDirections_ >
  inline constexpr bool isDual3rd< Dual3rd< nDirections_ > > = true;

  /** @brief Third-order vector-mode number: value, gradient, Hessian and third derivatives with respect to
   *  nDirections_ variables
   *  @tparam nDirections_ number of independent variables (directions)
   *
   *  The gradient and the Hessian are stored like in Dual2nd, the Hessian as its upper triangle, row by row. The third
   *  derivatives are stored for the sorted triples i <= j <= k of directions, ordered by i, then j, then k: the
   *  entries of a pair ( i, j ) form a row, which is contiguous in k.
   *
   *  The mask holds the directions the number depends on; the derivatives with respect to all other directions are
   *  zero.
   */
  template < std::size_t nDirections_ >
  class Dual3rd {
  public:
    /** @brief The second-order number with the same directions, whose layout of the gradient and the Hessian and whose
     *  loops over them this number uses
     */
    using Lower = Dual2nd< nDirections_ >;
    /** @brief Number of directions */
    static constexpr std::size_t nDirections = nDirections_;
    /** @brief Number of stored Hessian entries */
    static constexpr std::size_t nHessian = Lower::nHessian;
    /** @brief Number of stored third derivatives */
    static constexpr std::size_t nThird = nDirections * ( nDirections + 1 ) * ( nDirections + 2 ) / 6;
    static_assert( nDirections >= 1 && nDirections <= 64, "the direction mask has 64 bits" );
    /** @brief Mask of all directions */
    static constexpr std::uint64_t allDirections = Lower::allDirections;

    /** @brief Position of the third derivative ( i, j, j ), where the row of the pair ( i, j ), i <= j, starts
     *  @param i first direction
     *  @param j second direction, not less than i
     *  @return index into thirdDerivatives(); the entry ( i, j, k ), k >= j, is at this index plus k - j
     */
    static constexpr std::size_t thirdRowStart( std::size_t i, std::size_t j )
    {
      // the triples whose first direction is less than i, and the rows ( i, i ), ..., ( i, j - 1 )
      const std::size_t remaining = nDirections - i;
      return nThird - remaining * ( remaining + 1 ) * ( remaining + 2 ) / 6 + ( j - i ) * nDirections -
             ( j * ( j - 1 ) - i * ( i - 1 ) ) / 2;
    }

    /** @brief Position of the third derivative ( i, j, k ) among the stored sorted triples
     *  @param i first direction
     *  @param j second direction
     *  @param k third direction
     *  @return index into thirdDerivatives()
     */
    static constexpr std::size_t thirdIndex( std::size_t i, std::size_t j, std::size_t k )
    {
      if ( i > j )
        std::swap( i, j );
      if ( j > k )
        std::swap( j, k );
      if ( i > j )
        std::swap( i, j );
      return thirdRowStart( i, j ) + ( k - j );
    }

    /** @brief Calls f( i, j, q ) for every row of third derivatives, i.e., the entries ( i, j, k ), k >= j, of all
     *  pairs i <= j of directions in a mask; the rows of the other pairs do not change in the operations
     *  @tparam F callable
     *  @param mask directions
     *  @param f function of the two directions and the index q of the entry ( i, j, j ); the entry ( i, j, k ) is at
     *  q + k - j
     */
    template < typename F >
    [[gnu::always_inline]] static inline void forEachThirdRow( std::uint64_t mask, F&& f )
    {
      if ( mask == allDirections ) {
        std::size_t q = 0;
        for ( std::size_t i = 0; i < nDirections; i++ )
          for ( std::size_t j = i; j < nDirections; q += nDirections - j, j++ )
            f( i, j, q );
      }
      else
        for ( std::uint64_t iMask = mask; iMask; iMask &= iMask - 1 ) {
          const std::size_t i = static_cast< std::size_t >( std::countr_zero( iMask ) );
          for ( std::uint64_t jMask = iMask; jMask; jMask &= jMask - 1 ) {
            const std::size_t j = static_cast< std::size_t >( std::countr_zero( jMask ) );
            f( i, j, thirdRowStart( i, j ) );
          }
        }
    }

    /** @brief Calls f( q ) for every stored third derivative of the rows of a mask, see forEachThirdRow
     *  @tparam F callable
     *  @param mask directions
     *  @param f function of the index into thirdDerivatives()
     */
    template < typename F >
    [[gnu::always_inline]] static inline void forEachThirdEntry( std::uint64_t mask, F&& f )
    {
      if ( mask == allDirections )
        for ( std::size_t q = 0; q < nThird; q++ )
          f( q );
      else
        forEachThirdRow( mask, [&]( std::size_t, std::size_t j, std::size_t q ) {
          for ( std::size_t k = j; k < nDirections; k++, q++ )
            f( q );
        } );
    }

    /** @brief Tag of the constructor that leaves the derivatives uninitialized */
    struct Uninitialized {};

    /** @brief Creates a zero */
    Dual3rd() : g_{}, h_{}, t_{}, val_( 0.0 ), mask_( 0 ) {}

    /** @brief Creates a number with uninitialized derivatives, for the implementation of the operations, which must
     *  write all entries
     */
    explicit Dual3rd( Uninitialized ) : val_( 0.0 ), mask_( allDirections ) {}

    /** @brief Creates the result of an operation that depends on the directions in a mask: uninitialized if all
     *  directions, since the operation writes all entries then, and zero otherwise, since it writes only the
     *  entries of the rows in the mask
     *  @param mask directions of the result
     *  @return the number, with the mask set
     */
    [[gnu::always_inline]] static inline Dual3rd resultFor( std::uint64_t mask )
    {
      if ( mask == allDirections )
        return Dual3rd( Uninitialized{} );
      Dual3rd r;
      r.mask_ = mask;
      return r;
    }

    /** @brief Creates a constant, with zero derivatives
     *  @tparam S arithmetic type
     *  @param value value
     */
    template < ScalarMode::Constant S >
    Dual3rd( S value ) : g_{}, h_{}, t_{}, val_( static_cast< double >( value ) ), mask_( 0 )
    {
    }

    /** @brief Makes this number the independent variable `direction` (gradient = unit vector, zero higher
     *  derivatives)
     *  @param direction index of the variable
     */
    void seed( std::size_t direction )
    {
      if ( mask_ != 0 ) {
        g_.fill( 0.0 );
        h_.fill( 0.0 );
        t_.fill( 0.0 );
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
    double dxx( std::size_t i, std::size_t j ) const { return h_[Lower::hessianIndex( i, j )]; }

    /** @brief Third partial derivative
     *  @param i first direction
     *  @param j second direction
     *  @param k third direction
     *  @return third partial derivative with respect to variables i, j and k
     */
    double dxxx( std::size_t i, std::size_t j, std::size_t k ) const { return t_[thirdIndex( i, j, k )]; }

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

    /** @brief Third derivatives
     *  @return the third partial derivatives of the sorted triples of directions, see thirdIndex
     */
    const std::array< double, nThird >& thirdDerivatives() const { return t_; }

    /** @brief Third derivatives, for the implementation of the operations, which keep the mask consistent
     *  @return reference to the third partial derivatives of the sorted triples of directions, see thirdIndex
     */
    std::array< double, nThird >& thirdDerivatives() { return t_; }

    /** @brief Primal value
     *  @return value
     */
    explicit operator double() const { return val_; }

    /** @brief Adds a number
     *  @param b summand
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual3rd& operator+=( const Dual3rd& b )
    {
      combine( b, 1.0 );
      return *this;
    }

    /** @brief Subtracts a number
     *  @param b subtrahend
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual3rd& operator-=( const Dual3rd& b )
    {
      combine( b, -1.0 );
      return *this;
    }

    /** @brief Multiplies by a number
     *  @param b factor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual3rd& operator*=( const Dual3rd& b );

    /** @brief Divides by a number
     *  @param b divisor
     *  @return this number
     */
    [[gnu::always_inline]] inline Dual3rd& operator/=( const Dual3rd& b );

    /** @brief Adds a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual3rd& operator+=( S c )
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
    [[gnu::always_inline]] inline Dual3rd& operator-=( S c )
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
    [[gnu::always_inline]] inline Dual3rd& operator*=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ *= c_;
      Lower::forEachDirection( mask_, [&]( std::size_t i ) { g_[i] *= c_; } );
      Lower::forEachRow( mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] *= c_; } );
      forEachThirdEntry( mask_, [&]( std::size_t q ) { t_[q] *= c_; } );
      return *this;
    }

    /** @brief Divides by a constant
     *  @tparam S arithmetic type
     *  @param c constant
     *  @return this number
     */
    template < ScalarMode::Constant S >
    [[gnu::always_inline]] inline Dual3rd& operator/=( S c )
    {
      const double c_ = static_cast< double >( c );
      val_ /= c_;
      Lower::forEachDirection( mask_, [&]( std::size_t i ) { g_[i] /= c_; } );
      Lower::forEachRow( mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] /= c_; } );
      forEachThirdEntry( mask_, [&]( std::size_t q ) { t_[q] /= c_; } );
      return *this;
    }

  private:
    /** @brief this += sign b
     *  @param b number
     *  @param sign +1 or -1
     */
    [[gnu::always_inline]] inline void combine( const Dual3rd& b, double sign )
    {
      val_ += sign * b.val_;
      Lower::forEachDirection( b.mask_, [&]( std::size_t i ) { g_[i] += sign * b.g_[i]; } );
      Lower::forEachRow( b.mask_, [&]( std::size_t, std::size_t, std::size_t k ) { h_[k] += sign * b.h_[k]; } );
      forEachThirdEntry( b.mask_, [&]( std::size_t q ) { t_[q] += sign * b.t_[q]; } );
      mask_ |= b.mask_;
    }

    /** @brief Gradient */
    std::array< double, nDirections > g_;
    /** @brief Upper triangle of the Hessian, row by row */
    std::array< double, nHessian > h_;
    /** @brief Third derivatives of the sorted triples of directions */
    std::array< double, nThird > t_;
    /** @brief Value */
    double val_;
    /** @brief Directions this number depends on */
    std::uint64_t mask_;
  };

  /** @brief Applies a function, given by its value and first, second and third derivative at a.val(), to a number
   *  @tparam N number of directions
   *  @param a argument
   *  @param f0 f( a )
   *  @param f1 f'( a )
   *  @param f2 f''( a )
   *  @param f3 f'''( a )
   *  @return f( a ): gradient f1 grad a, Hessian f1 hess a + f2 grad a grad a^T, third derivatives
   *  f1 a_ijk + f2 ( a_i a_jk + a_j a_ik + a_k a_ij ) + f3 a_i a_j a_k
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > chain( const Dual3rd< N >& a, double f0, double f1, double f2, double f3 )
  {
    using D        = Dual3rd< N >;
    using L        = Dual2nd< N >;
    D r            = D::resultFor( a.mask() );
    r.val()        = f0;
    const auto& ag = a.gradient();
    const auto& ah = a.hessian();
    const auto& at = a.thirdDerivatives();
    auto&       rg = r.gradient();
    auto&       rh = r.hessian();
    auto&       rt = r.thirdDerivatives();
    L::forEachDirection( a.mask(), [&]( std::size_t i ) { rg[i] = f1 * ag[i]; } );
    L::forEachRow( a.mask(),
                   [&]( std::size_t i, std::size_t j, std::size_t k ) { rh[k] = f1 * ah[k] + f2 * ag[i] * ag[j]; } );
    D::forEachThirdRow( a.mask(), [&]( std::size_t i, std::size_t j, std::size_t q ) {
      // the Hessian entries ( i, k ) and ( j, k ), k >= j >= i, are at iRow + k and jRow + k
      const std::size_t iRow = L::rowStart( i ) - i, jRow = L::rowStart( j ) - j;
      const double      agi = ag[i], agj = ag[j], ahij = ah[iRow + j], f3gg = f3 * agi * agj;
      for ( std::size_t k = j; k < N; k++, q++ )
        rt[q] = f1 * at[q] + f2 * ( agi * ah[jRow + k] + agj * ah[iRow + k] + ag[k] * ahij ) + f3gg * ag[k];
    } );
    return r;
  }

  /** @brief Product: grad = a grad b + b grad a, hess = a hess b + b hess a + grad a grad b^T + grad b grad a^T, third
   *  derivatives a b_ijk + b a_ijk + a_i b_jk + a_j b_ik + a_k b_ij + b_i a_jk + b_j a_ik + b_k a_ij
   *  @tparam N number of directions
   *  @param a first factor
   *  @param b second factor
   *  @return a * b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > operator*( const Dual3rd< N >& a, const Dual3rd< N >& b )
  {
    using D                = Dual3rd< N >;
    using L                = Dual2nd< N >;
    const double        av = a.val(), bv = b.val();
    const std::uint64_t m = a.mask() | b.mask();
    D                   r = D::resultFor( m );
    r.val()               = av * bv;
    const auto& ag        = a.gradient();
    const auto& bg        = b.gradient();
    const auto& ah        = a.hessian();
    const auto& bh        = b.hessian();
    const auto& at        = a.thirdDerivatives();
    const auto& bt        = b.thirdDerivatives();
    auto&       rg        = r.gradient();
    auto&       rh        = r.hessian();
    auto&       rt        = r.thirdDerivatives();
    // rows outside m stay zero; within a row, the entries of the directions outside m compute to zero
    L::forEachDirection( m, [&]( std::size_t i ) { rg[i] = av * bg[i] + bv * ag[i]; } );
    L::forEachRow( m, [&]( std::size_t i, std::size_t j, std::size_t k ) {
      rh[k] = av * bh[k] + bv * ah[k] + ag[i] * bg[j] + bg[i] * ag[j];
    } );
    D::forEachThirdRow( m, [&]( std::size_t i, std::size_t j, std::size_t q ) {
      // the Hessian entries ( i, k ) and ( j, k ), k >= j >= i, are at iRow + k and jRow + k
      const std::size_t iRow = L::rowStart( i ) - i, jRow = L::rowStart( j ) - j;
      const double      agi = ag[i], agj = ag[j], ahij = ah[iRow + j];
      const double      bgi = bg[i], bgj = bg[j], bhij = bh[iRow + j];
      for ( std::size_t k = j; k < N; k++, q++ )
        rt[q] = av * bt[q] + bv * at[q] + agi * bh[jRow + k] + agj * bh[iRow + k] + ag[k] * bhij + bgi * ah[jRow + k] +
                bgj * ah[iRow + k] + bg[k] * ahij;
    } );
    return r;
  }

  template < std::size_t nDirections_ >
  Dual3rd< nDirections_ >& Dual3rd< nDirections_ >::operator*=( const Dual3rd& b )
  {
    return *this = *this * b;
  }

  /** @brief Reciprocal
   *  @tparam N number of directions
   *  @param a argument
   *  @return 1 / a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > reciprocal( const Dual3rd< N >& a )
  {
    const double r = 1.0 / a.val(), r2 = r * r;
    return chain( a, r, -r2, 2.0 * r2 * r, -6.0 * r2 * r2 );
  }

  /** @brief Quotient
   *  @tparam N number of directions
   *  @param a dividend
   *  @param b divisor
   *  @return a / b
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > operator/( const Dual3rd< N >& a, const Dual3rd< N >& b )
  {
    return a * reciprocal( b );
  }

  template < std::size_t nDirections_ >
  Dual3rd< nDirections_ >& Dual3rd< nDirections_ >::operator/=( const Dual3rd& b )
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
  [[gnu::always_inline]] inline Dual3rd< N > linearCombination( const Dual3rd< N >& a,
                                                                double              ca,
                                                                const Dual3rd< N >& b,
                                                                double              cb )
  {
    using D               = Dual3rd< N >;
    using L               = Dual2nd< N >;
    const std::uint64_t m = a.mask() | b.mask();
    D                   r = D::resultFor( m );
    r.val()               = ca * a.val() + cb * b.val();
    const auto& ag        = a.gradient();
    const auto& bg        = b.gradient();
    const auto& ah        = a.hessian();
    const auto& bh        = b.hessian();
    const auto& at        = a.thirdDerivatives();
    const auto& bt        = b.thirdDerivatives();
    auto&       rg        = r.gradient();
    auto&       rh        = r.hessian();
    auto&       rt        = r.thirdDerivatives();
    L::forEachDirection( m, [&]( std::size_t i ) { rg[i] = ca * ag[i] + cb * bg[i]; } );
    L::forEachRow( m, [&]( std::size_t, std::size_t, std::size_t k ) { rh[k] = ca * ah[k] + cb * bh[k]; } );
    D::forEachThirdEntry( m, [&]( std::size_t q ) { rt[q] = ca * at[q] + cb * bt[q]; } );
    return r;
  }

  /** @brief Product with a constant, writing every entry of the result once
   *  @tparam N number of directions
   *  @param a number
   *  @param c constant
   *  @return c a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > scaled( const Dual3rd< N >& a, double c )
  {
    using D        = Dual3rd< N >;
    using L        = Dual2nd< N >;
    D r            = D::resultFor( a.mask() );
    r.val()        = c * a.val();
    const auto& ag = a.gradient();
    const auto& ah = a.hessian();
    const auto& at = a.thirdDerivatives();
    auto&       rg = r.gradient();
    auto&       rh = r.hessian();
    auto&       rt = r.thirdDerivatives();
    L::forEachDirection( a.mask(), [&]( std::size_t i ) { rg[i] = c * ag[i]; } );
    L::forEachRow( a.mask(), [&]( std::size_t, std::size_t, std::size_t k ) { rh[k] = c * ah[k]; } );
    D::forEachThirdEntry( a.mask(), [&]( std::size_t q ) { rt[q] = c * at[q]; } );
    return r;
  }

  /** @brief Unary plus
   *  @tparam N number of directions
   *  @param a operand
   *  @return a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > operator+( const Dual3rd< N >& a )
  {
    return a;
  }

  /** @brief Negation
   *  @tparam N number of directions
   *  @param a operand
   *  @return -a
   */
  template < std::size_t N >
  [[gnu::always_inline]] inline Dual3rd< N > operator-( const Dual3rd< N >& a )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator+( const Dual3rd< N >& a, const Dual3rd< N >& b )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator-( const Dual3rd< N >& a, const Dual3rd< N >& b )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator+( Dual3rd< N > a, S c )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator+( S c, Dual3rd< N > a )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator-( Dual3rd< N > a, S c )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator-( S c, const Dual3rd< N >& a )
  {
    Dual3rd< N > r = scaled( a, -1.0 );
    return r += c;
  }

  /** @brief Product with a constant
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a number
   *  @param c constant
   *  @return a * c
   */
  template < std::size_t N, ScalarMode::Constant S >
  [[gnu::always_inline]] inline Dual3rd< N > operator*( const Dual3rd< N >& a, S c )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator*( S c, const Dual3rd< N >& a )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator/( Dual3rd< N > a, S c )
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
  [[gnu::always_inline]] inline Dual3rd< N > operator/( S c, const Dual3rd< N >& a )
  {
    const double c_ = static_cast< double >( c ), r = 1.0 / a.val(), r2 = r * r;
    return chain( a, c_ * r, -c_ * r2, 2.0 * c_ * r2 * r, -6.0 * c_ * r2 * r2 );
  }

  /** @brief An operand pair of a comparison involving at least one third-order vector-mode number
   *  @tparam A first operand type
   *  @tparam B second operand type
   */
  template < typename A, typename B >
  concept Dual3rdComparable = ( isDual3rd< A > && (isDual3rd< B > || ScalarMode::Constant< B >)) ||
                              ( ScalarMode::Constant< A > && isDual3rd< B > );

  /** @brief Compares the values
   *  @tparam A first operand type
   *  @tparam B second operand type
   *  @param a first operand
   *  @param b second operand
   *  @return value of a == value of b
   */
  template < typename A, typename B >
    requires Dual3rdComparable< A, B >
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
    requires Dual3rdComparable< A, B >
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
    requires Dual3rdComparable< A, B >
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
    requires Dual3rdComparable< A, B >
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
    requires Dual3rdComparable< A, B >
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
    requires Dual3rdComparable< A, B >
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
  Dual3rd< N > sqrt( const Dual3rd< N >& a )
  {
    const double x = a.val(), s = std::sqrt( x ), d = 0.5 / s;
    return chain( a, s, d, -0.5 * d / x, 0.75 * d / ( x * x ) );
  }

  /** @brief Cube root
   *  @tparam N number of directions
   *  @param a argument
   *  @return cbrt( a )
   */
  template < std::size_t N >
  Dual3rd< N > cbrt( const Dual3rd< N >& a )
  {
    const double x = a.val(), c = std::cbrt( x ), d = 1.0 / ( 3.0 * c * c );
    return chain( a, c, d, -2.0 / 3.0 * d / x, 10.0 / 9.0 * d / ( x * x ) );
  }

  /** @brief Exponential function
   *  @tparam N number of directions
   *  @param a argument
   *  @return exp( a )
   */
  template < std::size_t N >
  Dual3rd< N > exp( const Dual3rd< N >& a )
  {
    const double e = std::exp( a.val() );
    return chain( a, e, e, e, e );
  }

  /** @brief Natural logarithm
   *  @tparam N number of directions
   *  @param a argument
   *  @return log( a )
   */
  template < std::size_t N >
  Dual3rd< N > log( const Dual3rd< N >& a )
  {
    const double r = 1.0 / a.val();
    return chain( a, std::log( a.val() ), r, -r * r, 2.0 * r * r * r );
  }

  /** @brief Decimal logarithm
   *  @tparam N number of directions
   *  @param a argument
   *  @return log10( a )
   */
  template < std::size_t N >
  Dual3rd< N > log10( const Dual3rd< N >& a )
  {
    const double x = a.val(), r = 1.0 / ( x * std::numbers::ln10 );
    return chain( a, std::log10( x ), r, -r / x, 2.0 * r / ( x * x ) );
  }

  /** @brief Power with a constant exponent
   *  @tparam N number of directions
   *  @tparam S arithmetic type
   *  @param a base
   *  @param n exponent
   *  @return pow( a, n )
   */
  template < std::size_t N, ScalarMode::Constant S >
  Dual3rd< N > pow( const Dual3rd< N >& a, S n )
  {
    const double n_ = static_cast< double >( n ), x = a.val();
    if ( n_ == 2.0 )
      return chain( a, x * x, 2.0 * x, 2.0, 0.0 );
    if ( n_ == 0.0 )
      return Dual3rd< N >( 1.0 );
    if ( x != 0.0 ) {
      const double p = std::pow( x, n_ ), d = n_ * p / x, d2 = ( n_ - 1.0 ) * d / x;
      return chain( a, p, d, d2, ( n_ - 2.0 ) * d2 / x );
    }
    // a vanishing base: the derivatives without dividing by the base; those of an order above a positive integer
    // exponent vanish
    const double c2 = n_ * ( n_ - 1.0 ), c3 = c2 * ( n_ - 2.0 );
    return chain( a,
                  std::pow( x, n_ ),
                  n_ * std::pow( x, n_ - 1.0 ),
                  c2 == 0.0 ? 0.0 : c2 * std::pow( x, n_ - 2.0 ),
                  c3 == 0.0 ? 0.0 : c3 * std::pow( x, n_ - 3.0 ) );
  }

  /** @brief Power of a constant base
   *  @tparam S arithmetic type
   *  @tparam N number of directions
   *  @param b base
   *  @param a exponent
   *  @return pow( b, a )
   */
  template < ScalarMode::Constant S, std::size_t N >
  Dual3rd< N > pow( S b, const Dual3rd< N >& a )
  {
    const double b_ = static_cast< double >( b ), p = std::pow( b_, a.val() );
    const double l = b_ == 0.0 ? 0.0 : std::log( b_ );
    return chain( a, p, p * l, p * l * l, p * l * l * l );
  }

  /** @brief Power of two numbers
   *  @tparam N number of directions
   *  @param a base
   *  @param b exponent
   *  @return pow( a, b ) = exp( b log( a ) )
   */
  template < std::size_t N >
  Dual3rd< N > pow( const Dual3rd< N >& a, const Dual3rd< N >& b )
  {
    return exp( b * log( a ) );
  }

  /** @brief Sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return sin( a )
   */
  template < std::size_t N >
  Dual3rd< N > sin( const Dual3rd< N >& a )
  {
    const double s = std::sin( a.val() ), c = std::cos( a.val() );
    return chain( a, s, c, -s, -c );
  }

  /** @brief Cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return cos( a )
   */
  template < std::size_t N >
  Dual3rd< N > cos( const Dual3rd< N >& a )
  {
    const double s = std::sin( a.val() ), c = std::cos( a.val() );
    return chain( a, c, -s, -c, s );
  }

  /** @brief Tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return tan( a )
   */
  template < std::size_t N >
  Dual3rd< N > tan( const Dual3rd< N >& a )
  {
    const double t = std::tan( a.val() ), d = 1.0 + t * t;
    return chain( a, t, d, 2.0 * t * d, 2.0 * d * ( 1.0 + 3.0 * t * t ) );
  }

  /** @brief Arc sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return asin( a )
   */
  template < std::size_t N >
  Dual3rd< N > asin( const Dual3rd< N >& a )
  {
    const double x = a.val(), r = 1.0 / ( 1.0 - x * x ), d = std::sqrt( r );
    return chain( a, std::asin( x ), d, x * d * r, d * r * ( 1.0 + 3.0 * x * x * r ) );
  }

  /** @brief Arc cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return acos( a )
   */
  template < std::size_t N >
  Dual3rd< N > acos( const Dual3rd< N >& a )
  {
    const double x = a.val(), r = 1.0 / ( 1.0 - x * x ), d = std::sqrt( r );
    return chain( a, std::acos( x ), -d, -x * d * r, -d * r * ( 1.0 + 3.0 * x * x * r ) );
  }

  /** @brief Arc tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return atan( a )
   */
  template < std::size_t N >
  Dual3rd< N > atan( const Dual3rd< N >& a )
  {
    const double x = a.val(), d = 1.0 / ( 1.0 + x * x );
    return chain( a, std::atan( x ), d, -2.0 * x * d * d, 2.0 * d * d * ( 4.0 * x * x * d - 1.0 ) );
  }

  /** @brief Arc tangent of y / x, using the signs of both to determine the quadrant
   *  @tparam N number of directions
   *  @param y numerator
   *  @param x denominator
   *  @return atan2( y, x )
   */
  template < std::size_t N >
  Dual3rd< N > atan2( const Dual3rd< N >& y, const Dual3rd< N >& x )
  {
    // atan2 differs from atan( y / x ) and from -atan( x / y ) by piecewise constants only; the quotient with the
    // larger denominator is used
    Dual3rd< N > r = std::abs( x.val() ) >= std::abs( y.val() ) ? atan( y / x ) : -atan( x / y );
    r.val()        = std::atan2( y.val(), x.val() );
    return r;
  }

  /** @brief Hyperbolic sine
   *  @tparam N number of directions
   *  @param a argument
   *  @return sinh( a )
   */
  template < std::size_t N >
  Dual3rd< N > sinh( const Dual3rd< N >& a )
  {
    const double s = std::sinh( a.val() ), c = std::cosh( a.val() );
    return chain( a, s, c, s, c );
  }

  /** @brief Hyperbolic cosine
   *  @tparam N number of directions
   *  @param a argument
   *  @return cosh( a )
   */
  template < std::size_t N >
  Dual3rd< N > cosh( const Dual3rd< N >& a )
  {
    const double s = std::sinh( a.val() ), c = std::cosh( a.val() );
    return chain( a, c, s, c, s );
  }

  /** @brief Hyperbolic tangent
   *  @tparam N number of directions
   *  @param a argument
   *  @return tanh( a )
   */
  template < std::size_t N >
  Dual3rd< N > tanh( const Dual3rd< N >& a )
  {
    const double t = std::tanh( a.val() ), d = 1.0 - t * t;
    return chain( a, t, d, -2.0 * t * d, 2.0 * d * ( 3.0 * t * t - 1.0 ) );
  }

  /** @brief Absolute value; its derivatives at 0 are those of a
   *  @tparam N number of directions
   *  @param a argument
   *  @return |a|
   */
  template < std::size_t N >
  Dual3rd< N > abs( const Dual3rd< N >& a )
  {
    return a.val() < 0.0 ? -a : a;
  }

  /** @brief Absolute value; its derivatives at 0 are those of a
   *  @tparam N number of directions
   *  @param a argument
   *  @return |a|
   */
  template < std::size_t N >
  Dual3rd< N > fabs( const Dual3rd< N >& a )
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
  Dual3rd< N > max( const Dual3rd< N >& a, const Dual3rd< N >& b )
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
  Dual3rd< N > min( const Dual3rd< N >& a, const Dual3rd< N >& b )
  {
    return b.val() < a.val() ? b : a;
  }

  /** @brief True if the value is not a number
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isnan of the value
   */
  template < std::size_t N >
  bool isnan( const Dual3rd< N >& a )
  {
    return std::isnan( a.val() );
  }

  /** @brief True if the value is infinite
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isinf of the value
   */
  template < std::size_t N >
  bool isinf( const Dual3rd< N >& a )
  {
    return std::isinf( a.val() );
  }

  /** @brief True if the value is finite
   *  @tparam N number of directions
   *  @param a number
   *  @return std::isfinite of the value
   */
  template < std::size_t N >
  bool isfinite( const Dual3rd< N >& a )
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
  std::ostream& operator<<( std::ostream& out, const Dual3rd< N >& a )
  {
    return out << a.val();
  }

} // namespace Marmot::Math::Differentiation::Automatic::VectorMode

namespace Eigen {
  /** @brief Eigen support for third-order vector-mode numbers
   *  @tparam N number of directions
   */
  template < std::size_t N >
  struct NumTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N > > : NumTraits< double > {
    /** @brief Real type */
    using Real = Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N >;
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

  /** @brief Lets Eigen mix third-order vector-mode numbers with double
   *  @tparam N number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t N, typename BinaryOp >
  struct ScalarBinaryOpTraits< Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N >, double, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N >;
  };

  /** @brief Lets Eigen mix double with third-order vector-mode numbers
   *  @tparam N number of directions
   *  @tparam BinaryOp Eigen binary functor
   */
  template < std::size_t N, typename BinaryOp >
  struct ScalarBinaryOpTraits< double, Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N >, BinaryOp > {
    /** @brief Result type */
    using ReturnType = Marmot::Math::Differentiation::Automatic::VectorMode::Dual3rd< N >;
  };
} // namespace Eigen
