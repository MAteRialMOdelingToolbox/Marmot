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
 *  @brief Backend-neutral number types and primitive operations for forward-mode automatic differentiation.
 *
 *  This is the only Marmot header that names the AD backends. All other Marmot code uses the types and functions
 *  declared here, so that the backends can be exchanged in one place.
 *
 *  Backends:
 *  - autodiff (https://autodiff.github.io), always used: HyperDual, Taylor and, by default, JacobianDual.
 *  - Sacado (https://trilinos.github.io/sacado.html), optional: if Marmot is configured with
 *    MARMOT_AD_VECTOR_BACKEND=Sacado (which defines MARMOT_AD_VECTOR_BACKEND_SACADO), JacobianDual is a
 *    Sacado::Fad::SFad, and the first-derivative drivers (df_dT, dF_dT, fixed-size dF_dX) compute all directions in a
 *    single evaluation (vector mode) instead of one direction per evaluation.
 *
 *  Vocabulary:
 *  - HyperDual<order>: nested dual number with one independent direction per nesting level; level 1 is the
 *    outermost level. dual, dual2nd and dual3rd are the orders 1, 2 and 3.
 *  - Taylor<order>: univariate truncated Taylor number carrying the derivatives 0..order of one variable.
 *  - JacobianDual<nDirections>: number type of the first-derivative drivers df_dT, dF_dT and dF_dX< n >, with
 *    nDirections independent variables. Code differentiated by these drivers must use this type, not dual, to work
 *    with both backends. With autodiff it is dual.
 *  - seed<level>( x, v ) sets the derivative of x along nesting level `level` to v, derivative<level>( y ) reads the
 *    mixed derivative of y along the levels 1..level, value( y ) reads the primal value.
 *  - For vector-mode numbers (isVectorDual), seedDirection( x, i ) makes x the independent variable i and
 *    partialDerivative( y, i ) reads the partial derivative of y with respect to variable i.
 */

#pragma once
#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <cstddef>
#include <type_traits>
#ifdef MARMOT_AD_VECTOR_BACKEND_SACADO
#  include <Sacado_No_Kokkos.hpp>
#endif

namespace Marmot::Math::Differentiation::Automatic {

  /** @brief Hyper-dual number of the given order
   *  @tparam order number of nesting levels, i.e., the highest derivative order it carries
   *
   *  A hyper-dual number of order N holds 2^N coefficients: the primal value, the first derivatives along each of the
   *  N nesting levels and all mixed derivatives. Seeding different variables on different levels yields mixed
   *  higher-order derivatives in one evaluation.
   */
  template < size_t order >
  using HyperDual = autodiff::HigherOrderDual< order, double >;

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
  using Taylor = autodiff::Real< order, double >;

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
  inline constexpr bool isHyperDual = autodiff::detail::isDual< Scalar >;

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
    autodiff::seed< level >( x, static_cast< double >( value ) );
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
    return autodiff::derivative< level >( x );
  }

  /** @brief Reads the primal value of a hyper-dual or vector-mode number
   *  @tparam Scalar hyper-dual number type or JacobianDual
   *  @param x number
   *  @return primal (zeroth order) value of x
   */
  template < typename Scalar >
  double value( const Scalar& x );

  namespace Detail {
    /** @brief Implementation of isVectorDual
     *  @tparam Scalar type to check
     */
    template < typename Scalar >
    struct IsVectorDual : std::false_type {};
#ifdef MARMOT_AD_VECTOR_BACKEND_SACADO
    /** @brief Implementation of isVectorDual for Sacado::Fad::SFad
     *  @tparam nDirections number of derivatives
     */
    template < int nDirections >
    struct IsVectorDual< Sacado::Fad::SFad< double, nDirections > > : std::true_type {};
#endif
  } // namespace Detail

  /** @brief Number type of the first-derivative drivers (df_dT, dF_dT, dF_dX< n >)
   *  @tparam nDirections number of independent variables
   *
   *  With the Sacado backend a Sacado::Fad::SFad< double, nDirections >, which carries all partial derivatives after a
   *  single evaluation; otherwise dual, which carries one direction per evaluation.
   */
#ifdef MARMOT_AD_VECTOR_BACKEND_SACADO
  template < size_t nDirections >
  using JacobianDual = Sacado::Fad::SFad< double, static_cast< int >( nDirections ) >;
#else
  template < size_t nDirections >
  using JacobianDual = dual;
#endif

  /** @brief True if Scalar is a vector-mode number, i.e., a JacobianDual of the Sacado backend
   *  @tparam Scalar type to check
   */
  template < typename Scalar >
  inline constexpr bool isVectorDual = Detail::IsVectorDual< Scalar >::value;

  template < typename Scalar >
  double value( const Scalar& x )
  {
    if constexpr ( isVectorDual< Scalar > )
      return x.val();
    else
      return autodiff::detail::val( x );
  }

  /** @brief Makes a vector-mode number the independent variable `direction`
   *  @tparam Scalar vector-mode number type (isVectorDual)
   *  @param x number to seed; its derivatives must be zero
   *  @param direction index of the independent variable, less than the number of directions of Scalar
   */
  template < typename Scalar >
  void seedDirection( Scalar& x, size_t direction )
  {
    x.fastAccessDx( static_cast< int >( direction ) ) = 1.0;
  }

  /** @brief Reads a partial derivative of a vector-mode number
   *  @tparam Scalar vector-mode number type (isVectorDual)
   *  @param x number
   *  @param i index of the independent variable
   *  @return partial derivative of x with respect to variable i
   */
  template < typename Scalar >
  double partialDerivative( const Scalar& x, size_t i )
  {
    return x.fastAccessDx( static_cast< int >( i ) );
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
    Taylor< order >                 x_( x );
    const auto                      d = autodiff::derivatives( f, autodiff::along( 1. ), autodiff::at( x_ ) );
    std::array< double, order + 1 > out{};
    for ( size_t k = 0; k <= order; ++k )
      out[k] = d[k];
    return out;
  }

} // namespace Marmot::Math::Differentiation::Automatic

#ifdef MARMOT_AD_VECTOR_BACKEND_SACADO
namespace Eigen {
  /** @brief Eigen support for the vector-mode numbers of the Sacado backend */
  template < int nDirections >
  struct NumTraits< Sacado::Fad::SFad< double, nDirections > > : NumTraits< double > {
    using Real       = Sacado::Fad::SFad< double, nDirections >;
    using NonInteger = Real;
    using Literal    = double;
    using Nested     = Real;
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

  /** @brief Lets Eigen mix the vector-mode numbers of the Sacado backend with double */
  template < int nDirections, typename BinaryOp >
  struct ScalarBinaryOpTraits< Sacado::Fad::SFad< double, nDirections >, double, BinaryOp > {
    using ReturnType = Sacado::Fad::SFad< double, nDirections >;
  };

  /** @brief Lets Eigen mix double with the vector-mode numbers of the Sacado backend */
  template < int nDirections, typename BinaryOp >
  struct ScalarBinaryOpTraits< double, Sacado::Fad::SFad< double, nDirections >, BinaryOp > {
    using ReturnType = Sacado::Fad::SFad< double, nDirections >;
  };
} // namespace Eigen
#endif
