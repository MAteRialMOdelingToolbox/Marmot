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

/** @file MarmotAutomaticDifferentiationTypes.h
 *  @brief Number types and primitive operations for forward-mode automatic differentiation.
 *
 *  All Marmot code uses the types and functions declared here instead of naming the implementations directly:
 *  - Hyper-dual and Taylor numbers (dual, dual2nd, dual3rd, HyperDual, Taylor): ScalarMode
 *    (MarmotAutomaticDifferentiationScalarMode.h).
 *  - JacobianDual, the number type of the first-derivative drivers (df_dT, dF_dT, fixed-size dF_dX): VectorMode::Dual
 *    (MarmotAutomaticDifferentiationVectorMode.h), which computes all directions in a single evaluation.
 *  - HessianDual, the number type of the second-derivative driver d2f_dT2: VectorMode::Dual2nd
 *    (MarmotAutomaticDifferentiationVectorModeSecondOrder.h), which computes the gradient and the Hessian in a single
 *    evaluation.
 *
 *  Vocabulary:
 *  - HyperDual<order>: nested dual number with one independent direction per nesting level; level 1 is the
 *    outermost level. dual, dual2nd and dual3rd are the orders 1, 2 and 3.
 *  - Taylor<order>: univariate truncated Taylor number carrying the derivatives 0..order of one variable.
 *  - JacobianDual<nDirections>: number type of the first-derivative drivers df_dT, dF_dT and dF_dX< n >, with
 *    nDirections independent variables. Code differentiated by these drivers must use this type, not dual.
 *  - seed<level>( x, v ) sets the derivative of x along nesting level `level` to v, derivative<level>( y ) reads the
 *    mixed derivative of y along the levels 1..level, value( y ) reads the primal value.
 *  - HessianDual<nDirections>: number type of the second-derivative driver d2f_dT2, with nDirections independent
 *    variables. Code differentiated by d2f_dT2 must use this type, not dual2nd.
 *  - For vector-mode numbers (isVectorDual, isHessianDual), seedDirection( x, i ) makes x the independent variable i
 *    and partialDerivative( y, i ) reads the partial derivative of y with respect to variable i;
 *    secondPartialDerivative( y, i, j ) reads a second partial derivative of a HessianDual.
 */

#pragma once
#include "Marmot/MarmotAutomaticDifferentiationScalarMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorMode.h"
#include "Marmot/MarmotAutomaticDifferentiationVectorModeSecondOrder.h"
#include <array>
#include <cstddef>
#include <type_traits>

namespace Marmot::Math::Differentiation::Automatic {

  /** @brief Hyper-dual number of the given order
   *  @tparam order number of nesting levels, i.e., the highest derivative order it carries
   *
   *  A hyper-dual number of order N holds 2^N coefficients: the primal value, the first derivatives along each of the
   *  N nesting levels and all mixed derivatives. Seeding different variables on different levels yields mixed
   *  higher-order derivatives in one evaluation.
   */
  template < size_t order >
  using HyperDual = ScalarMode::HyperDual< order >;

  /** @brief First order dual number */
  using dual = HyperDual< 1 >;

  /** @brief Second order hyper-dual number */
  using dual2nd = HyperDual< 2 >;

  /** @brief Third order hyper-dual number */
  using dual3rd = HyperDual< 3 >;

  /** @brief Univariate truncated Taylor number carrying the derivatives 0..order of one variable
   *  @tparam order highest derivative order
   */
  template < size_t order >
  using Taylor = ScalarMode::Taylor< order >;

  /** @brief Dynamic-size Eigen column vector of first order dual numbers */
  using VectorXdual = Eigen::Matrix< dual, Eigen::Dynamic, 1 >;

  /** @brief Dynamic-size Eigen matrix of first order dual numbers */
  using MatrixXdual = Eigen::Matrix< dual, Eigen::Dynamic, Eigen::Dynamic >;

  /** @brief Dynamic-size Eigen column vector of second order hyper-dual numbers */
  using VectorXdual2nd = Eigen::Matrix< dual2nd, Eigen::Dynamic, 1 >;

  /** @brief True if Scalar is a hyper-dual number of any order
   *  @tparam Scalar type to check
   */
  template < typename Scalar >
  inline constexpr bool isHyperDual = ScalarMode::isDual< Scalar >;

  /** @brief Sets the derivative of a hyper-dual number along one nesting level
   *  @tparam level nesting level (1 = outermost) to seed
   *  @tparam Scalar hyper-dual number type
   *  @tparam Value type of the seed value
   *  @param x hyper-dual number to seed
   *  @param value seed value, usually 1.0 to make x the independent variable of this level, or 0.0 to unseed it
   *
   *  Seeding overwrites the whole derivative node of the level, i.e., also the mixed derivatives it contains.
   */
  template < size_t level, typename Scalar, typename Value >
  void seed( Scalar& x, Value&& value )
  {
    static_assert( level >= 1 && level <= ScalarMode::orderOf< Scalar >, "invalid nesting level" );
    if constexpr ( level == 1 )
      x.grad() = static_cast< double >( value );
    else
      seed< level - 1 >( x.val(), value );
  }

  /** @brief Reads the mixed derivative of a hyper-dual number along the nesting levels 1..level
   *  @tparam level highest nesting level of the mixed derivative
   *  @tparam Scalar hyper-dual number type
   *  @param x hyper-dual number
   *  @return derivative of x along the levels 1, ..., level
   */
  template < size_t level, typename Scalar >
  double derivative( const Scalar& x )
  {
    if constexpr ( level == 0 )
      return ScalarMode::primalValue( x );
    else
      return derivative< level - 1 >( x.grad() );
  }

  /** @brief Reads the primal value of a hyper-dual or vector-mode number
   *  @tparam Scalar hyper-dual number type or JacobianDual
   *  @param x number
   *  @return primal (zeroth order) value of x
   */
  template < typename Scalar >
  double value( const Scalar& x );

  /** @brief Number type of the first-derivative drivers (df_dT, dF_dT, dF_dX< n >)
   *  @tparam nDirections number of independent variables
   *
   *  A VectorMode::Dual< nDirections >, which carries all partial derivatives after a single evaluation.
   */
  template < size_t nDirections >
  using JacobianDual = VectorMode::Dual< nDirections >;

  /** @brief Number type of the second-derivative driver d2f_dT2
   *  @tparam nDirections number of independent variables
   *
   *  A VectorMode::Dual2nd< nDirections >, which carries the gradient and the Hessian after a single evaluation.
   */
  template < size_t nDirections >
  using HessianDual = VectorMode::Dual2nd< nDirections >;

  /** @brief True if Scalar is a second-order vector-mode number (VectorMode::Dual2nd)
   *  @tparam Scalar type to check
   */
  template < typename Scalar >
  inline constexpr bool isHessianDual = VectorMode::isDual2nd< Scalar >;

  /** @brief True if Scalar is a first-order vector-mode number (VectorMode::Dual)
   *  @tparam Scalar type to check
   */
  template < typename Scalar >
  inline constexpr bool isVectorDual = VectorMode::isDual< Scalar >;

  template < typename Scalar >
  double value( const Scalar& x )
  {
    if constexpr ( isVectorDual< Scalar > || isHessianDual< Scalar > )
      return x.val();
    else
      return ScalarMode::primalValue( x );
  }

  /** @brief Makes a vector-mode number the independent variable `direction`
   *  @tparam Scalar vector-mode number type (isVectorDual or isHessianDual)
   *  @param x number to seed; its derivatives must be zero
   *  @param direction index of the independent variable, less than the number of directions of Scalar
   */
  template < typename Scalar >
  void seedDirection( Scalar& x, size_t direction )
  {
    if constexpr ( isHessianDual< Scalar > )
      x.seed( direction );
    else
      x.dx( direction ) = 1.0;
  }

  /** @brief Reads a partial derivative of a vector-mode number
   *  @tparam Scalar vector-mode number type (isVectorDual or isHessianDual)
   *  @param x number
   *  @param i index of the independent variable
   *  @return partial derivative of x with respect to variable i
   */
  template < typename Scalar >
  double partialDerivative( const Scalar& x, size_t i )
  {
    return x.dx( i );
  }

  /** @brief Reads a second partial derivative of a second-order vector-mode number
   *  @tparam Scalar second-order vector-mode number type (isHessianDual)
   *  @param x number
   *  @param i index of the first independent variable
   *  @param j index of the second independent variable
   *  @return second partial derivative of x with respect to variables i and j
   */
  template < typename Scalar >
  double secondPartialDerivative( const Scalar& x, size_t i, size_t j )
  {
    return x.dxx( i, j );
  }

  /** @brief Contiguous coefficient storage of a hyper-dual number
   *  @tparam order order of the hyper-dual number
   *  @param x hyper-dual number
   *  @return pointer to the 2^order coefficients of x
   *
   *  The coefficients are stored level by level: the first half holds the hyper-dual number of order - 1 that is
   *  the value node of the outermost level, the second half the one that is its derivative node. Hence the primal value
   *  is at index 0 and the first derivative along level k at index 2^(order - k). This layout is what the
   * order-shifting functions in MarmotAutomaticDifferentiation.h rely on.
   */
  template < size_t order >
  double* coefficients( HyperDual< order >& x )
  {
    static_assert( std::is_standard_layout_v< HyperDual< order > > &&
                     sizeof( HyperDual< order > ) == ( size_t( 1 ) << order ) * sizeof( double ),
                   "HyperDual must be a dense array of 2^order doubles" );
    return reinterpret_cast< double* >( &x );
  }

  /** @brief Contiguous coefficient storage of a hyper-dual number (read-only)
   *  @tparam order order of the hyper-dual number
   *  @param x hyper-dual number
   *  @return pointer to the 2^order coefficients of x
   */
  template < size_t order >
  const double* coefficients( const HyperDual< order >& x )
  {
    return coefficients< order >( const_cast< HyperDual< order >& >( x ) );
  }

  /** @brief Computes the derivatives 0..order of a univariate function with Taylor-mode automatic differentiation
   *  @tparam order highest derivative order
   *  @tparam Function callable mapping a Taylor<order> to a Taylor<order>
   *  @param f function to differentiate
   *  @param x point at which the derivatives are evaluated
   *  @return the value and the derivatives of f at x; entry k holds the k-th derivative
   */
  template < size_t order, typename Function >
  std::array< double, order + 1 > taylorDerivatives( const Function& f, double x )
  {
    std::array< double, order + 1 > out{};
    Taylor< order >                 x_( x );
    if constexpr ( order >= 1 )
      x_[1] = 1.0;
    const Taylor< order > y         = f( x_ );
    double                factorial = 1.0;
    for ( size_t k = 0; k <= order; ++k ) {
      factorial *= k ? double( k ) : 1.0;
      out[k] = y[k] * factorial;
    }
    return out;
  }

} // namespace Marmot::Math::Differentiation::Automatic
