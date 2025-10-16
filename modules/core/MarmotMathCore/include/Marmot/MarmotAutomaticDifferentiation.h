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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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

#pragma once
#include "autodiff/forward/dual/eigen.hpp"
#include <autodiff/forward/dual/dual.hpp>
#include <functional>

namespace Marmot {

  namespace AutomaticDifferentiation {

    using namespace autodiff;
    using namespace Eigen;

    /** @brief Creates a 2nd order hyper-dual number from a dual number
     *  @param x Input dual number
     *  @return 2nd order hyper-dual number
     *
     *  This function takes a dual number as input and returns a 2nd order hyper-dual number.
     *  This is done by keeping the real part and shifting the first derivative part to the second derivative part.
     */
    dual2nd shiftTo2ndOrderDual( const dual& x );

    /** @brief Creates a vector of 2nd order hyper-dual numbers from a vector of dual numbers
     *  @param X Input vector of dual numbers
     *  @return Vector of 2nd order hyper-dual numbers
     *
     *  This function takes a vector of dual numbers as input and returns a vector of 2nd order hyper-dual numbers.
     *  This is done by keeping the real part and shifting the first derivative part to the second derivative part for
     * each element in the vector.
     */
    VectorXdual2nd shiftTo2ndOrderDual( const VectorXdual& X );

    /** @brief Access the value node of a dual number at a specific order
     *  @tparam order The order of the dual number to access
     *  @tparam T The type of the value stored in the dual number
     *  @tparam G The type of the gradient stored in the dual number
     *  @param dual The dual number to access
     *  @return A reference to the value node at the specified order
     *
     *  This function recursively accesses the value node of a dual number at a specific order.
     *  It uses template metaprogramming to ensure that the order is valid for the given dual number type.
     *  If the order is 0, it returns the value of the dual number. If the order is greater than 0,
     *  it recursively calls itself on the value of the dual number until it reaches the desired order.
     */
    template < size_t order, typename T, typename G >
    const auto& valnode( const Dual< T, G >& dual )
    {
      constexpr auto N = detail::Order< Dual< T, G > >;
      static_assert( order <= N );
      if constexpr ( order == 0 )
        return dual.val;
      else if constexpr ( order == 1 )
        return dual.val;
      else
        return valnode< order - 1 >( dual.val );
    }

    template < size_t order, typename T, typename G >
    auto& valnode( Dual< T, G >& dual )
    {
      constexpr auto N = detail::Order< Dual< T, G > >;
      static_assert( order <= N );
      if constexpr ( order == 0 )
        return dual.val;
      else if constexpr ( order == 1 )
        return dual.val;
      else
        return valnode< order - 1 >( dual.val );
    }

    /** @brief Increases the order of a hyper-dual number by 1
     *  @tparam order The current order of the dual number
     *  @param in The input dual number
     *  @return A new dual number with order increased by 1
     *
     *  This function takes a dual number of a given order and creates a new dual number with the order increased by 1.
     *  The derivative parts are elevated (shifted) one order.
     */
    template < size_t order >
    autodiff::HigherOrderDual< order + 1, double > increaseDualOrderWithShift(
      const autodiff::HigherOrderDual< order, double >& in )
    {
      using out_scalar_type = autodiff::HigherOrderDual< order + 1, double >;
      using namespace autodiff::detail;

      out_scalar_type out( 0.0 );
      const double*   in_point  = &valnode< order >( in );
      double*         out_point = &valnode< order + 1 >( out );

      for ( size_t i = 0; i < size_t( std::pow( 2, order ) ); i++ ) {
        *( out_point + i ) = *( in_point + i );
      }

      return out;
    }

    /** @brief Decreases the order of a hyper-dual number by 1
     *  @tparam order The current order of the dual number
     *  @param in The input dual number
     *  @return A new dual number with order decreased by 1
     *
     *  This function takes a dual number of a given order and creates a new dual number with the order decreased by 1.
     *  The derivative parts are copied to the lower order dual number.
     *  The highest derivative part of the input dual number is discarded.
     */
    template < size_t order >
    autodiff::HigherOrderDual< order - 1, double > decreaseDualOrder( autodiff::HigherOrderDual< order, double >& in )
    {
      using out_scalar_type = autodiff::HigherOrderDual< order - 1, double >;
      using namespace autodiff::detail;

      out_scalar_type out( 0.0 );
      double*         in_point  = &valnode< order >( in );
      double*         out_point = &valnode< order - 1 >( out );

      for ( size_t i = 0; i < size_t( std::pow( 2, order - 1 ) ); i++ ) {
        *( out_point + i ) = *( in_point + i );
      }

      return out;
    }

    /** @brief Decreases the order of a hyper-dual number by 1 with a shift
     *  @tparam order The current order of the dual number
     *  @param in The input dual number
     *  @return A new dual number with order decreased by 1
     *
     *  This function takes a dual number of a given order and creates a new dual number with the order decreased by 1.
     *  The derivative parts are shifted down by one order.
     *  The first derivative part of the input dual number is discarded.
     */
    template < size_t order >
    autodiff::HigherOrderDual< order - 1, double > decreaseDualOrderWithShift(
      autodiff::HigherOrderDual< order, double >& in )
    {
      using out_scalar_type = autodiff::HigherOrderDual< order - 1, double >;
      using namespace autodiff::detail;

      out_scalar_type out( 0.0 );
      double*         in_point  = &valnode< order >( in );
      double*         out_point = &valnode< order - 1 >( out );

      for ( size_t i = 0; i < size_t( std::pow( 2, order - 1 ) ); i++ ) {
        *( out_point + i ) = *( in_point + i + size_t( std::pow( 2, order - 1 ) ) );
      }

      return out;
    }

    /** @brief Increases the order of a vector of hyper-dual numbers by 1
     *  @tparam order The current order of the dual numbers in the vector
     *  @param in The input vector of dual numbers
     *  @return A new vector of dual numbers with order increased by 1
     *
     *  This function takes a vector of dual numbers of a given order and creates a new vector of dual numbers with the
     * order increased by 1. The derivative parts are elevated (shifted) one order for each element in the vector.
     */
    template < size_t order >
    Vector< HigherOrderDual< order + 1, double >, -1 > increaseDualOrderWithShift(
      const Vector< HigherOrderDual< order, double >, -1 >& in )
    {
      using in_scalar_type  = HigherOrderDual< order, double >;
      using out_scalar_type = HigherOrderDual< order + 1, double >;

      Vector< out_scalar_type, -1 > out      = Vector< out_scalar_type, -1 >( in.size() );
      out_scalar_type*              out_data = out.data();
      const in_scalar_type*         in_data  = in.data();

      for ( int i = 0; i < in.size(); i++ ) {
        out_data[i] = increaseDualOrderWithShift< order >( in_data[i] );
      }
      return out;
    }

    /** @typedef scalar_to_scalar_function_type
     *  @brief A type alias for a scalar-to-scalar function that takes and returns dual numbers
     */
    using scalar_to_scalar_function_type = std::function< dual( const dual& ) >;

    /** @brief Computes the derivative of a scalar-to-scalar function at a given point using automatic differentiation
     *  @param f The scalar-to-scalar function to differentiate
     *  @param x The point at which to evaluate the derivative
     *  @return The derivative of the function at the given point
     *
     *  This function takes a scalar-to-scalar function and a point as input, and returns the derivative of the function
     * at that point. It uses automatic differentiation to compute the derivative accurately.
     */
    double df_dx( const scalar_to_scalar_function_type& f, const double& x );

    /** @typedef scalar_to_scalar_function_type_2nd
     *  @brief A type alias for a scalar-to-scalar function that takes and returns second order dual numbers
     */
    using scalar_to_scalar_function_type_2nd = std::function< dual2nd( const dual2nd& ) >;

    /** @brief Computes the derivative of a dual-valued scalar-to-scalar function at a given point using automatic
     * differentiation
     *  @param f The scalar-to-scalar function to differentiate
     *  @param x The point at which to evaluate the derivative
     *  @return The derivative of the function at the given point
     *
     *  This function takes a 2nd order dual-valued scalar-to-scalar function and a point as input, and returns the
     * derivative of the function at that point. Note that the function is 2nd order dual-valued, but the input is only
     * 1st order dual-valued.
     */
    dual df_dx( const scalar_to_scalar_function_type_2nd& f, const dual& x );

    /** @typedef vector_to_vector_function_type_dual
     *  @brief A type alias for a vector-to-vector function that takes and returns dual-valued vectors
     */
    using vector_to_vector_function_type_dual = std::function< VectorXdual( const VectorXdual& X ) >;

    /** @brief Computes the Jacobian of a vector-to-vector function at a given point using automatic differentiation
     *  @param F The vector-to-vector function to differentiate
     *  @param X The point at which to evaluate the Jacobian
     *  @return A pair containing the function value and the Jacobian matrix at the given point
     *
     *  This function takes a vector-to-vector function and a point as input, and returns a pair containing the function
     * value and the Jacobian matrix at that point. It uses automatic differentiation to compute the Jacobian
     * accurately.
     */
    std::pair< VectorXd, MatrixXd > dF_dX( const vector_to_vector_function_type_dual& F, const VectorXd& X );

    /** @typedef vector_to_vector_function_type_dual2nd
     *  @brief A type alias for a vector-to-vector function that takes and returns second order dual-valued vectors
     */
    using vector_to_vector_function_type_dual2nd = std::function< VectorXdual2nd( const VectorXdual2nd& X ) >;

    /** @brief Computes the Jacobian of a dual-valued vector-to-vector function at a given point using automatic
     * differentiation
     *  @param F The vector-to-vector function to differentiate
     *  @param X The point at which to evaluate the Jacobian
     *  @return A pair containing the function value and the Jacobian matrix at the given point
     *
     *  This function takes a 2nd order dual-valued vector-to-vector function and a point as input, and returns a
     * dual-valued pair containing the function value and the Jacobian matrix at that point. Note that the function is
     * 2nd order dual-valued, but the input is only 1st order dual-valued.
     */
    std::pair< VectorXdual, MatrixXdual > dF_dX_2nd( const vector_to_vector_function_type_dual2nd& F,
                                                     const VectorXdual&                            X );
  } // namespace AutomaticDifferentiation

} // namespace Marmot
