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

#pragma once
#include "Marmot/MarmotAutomaticDifferentiationTypes.h"
#include <Eigen/Core>
#include <functional>

namespace Marmot::Math::Differentiation::Automatic {

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

  /** @brief Increases the order of a hyper-dual number by 1
   *  @tparam order The current order of the dual number
   *  @param in The input dual number
   *  @return A new dual number with order increased by 1
   *
   *  This function takes a dual number of a given order and creates a new dual number with the order increased by 1.
   *  The derivative parts are elevated (shifted) one order.
   */
  template < size_t order >
  HyperDual< order + 1 > increaseDualOrderWithShift( const HyperDual< order >& in )
  {
    HyperDual< order + 1 > out( 0.0 );
    const double*          in_point  = coefficients< order >( in );
    double*                out_point = coefficients< order + 1 >( out );

    for ( size_t i = 0; i < ( size_t( 1 ) << order ); i++ ) {
      out_point[i] = in_point[i];
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
  HyperDual< order - 1 > decreaseDualOrder( HyperDual< order >& in )
  {
    HyperDual< order - 1 > out( 0.0 );
    const double*          in_point  = coefficients< order >( in );
    double*                out_point = coefficients< order - 1 >( out );

    for ( size_t i = 0; i < ( size_t( 1 ) << ( order - 1 ) ); i++ ) {
      out_point[i] = in_point[i];
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
  HyperDual< order - 1 > decreaseDualOrderWithShift( HyperDual< order >& in )
  {
    constexpr size_t       half = size_t( 1 ) << ( order - 1 );
    HyperDual< order - 1 > out( 0.0 );
    const double*          in_point  = coefficients< order >( in );
    double*                out_point = coefficients< order - 1 >( out );

    for ( size_t i = 0; i < half; i++ ) {
      out_point[i] = in_point[i + half];
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
  Eigen::Vector< HyperDual< order + 1 >, -1 > increaseDualOrderWithShift(
    const Eigen::Vector< HyperDual< order >, -1 >& in )
  {
    using in_scalar_type  = HyperDual< order >;
    using out_scalar_type = HyperDual< order + 1 >;

    Eigen::Vector< out_scalar_type, -1 > out      = Eigen::Vector< out_scalar_type, -1 >( in.size() );
    out_scalar_type*                     out_data = out.data();
    const in_scalar_type*                in_data  = in.data();

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
  std::pair< Eigen::VectorXd, Eigen::MatrixXd > dF_dX( const vector_to_vector_function_type_dual& F,
                                                       const Eigen::VectorXd&                     X );

  /** @typedef fixed_size_vector_to_vector_function_type
   *  @brief A type alias for a function mapping a fixed-size vector to a vector of the same size, with JacobianDual
   *  numbers (one direction per vector component)
   *  @tparam n size of the vectors
   */
  template < int n >
  using fixed_size_vector_to_vector_function_type = std::function< Eigen::Vector< JacobianDual< n >, n >(
    const Eigen::Vector< JacobianDual< n >, n >& X ) >;

  /** @brief Computes the Jacobian of a fixed-size vector-to-vector function using automatic differentiation
   *  @tparam n size of the vectors; not deduced from a lambda, so it must be given explicitly, e.g. dF_dX< 6 >( ... )
   *  @param F The vector-to-vector function to differentiate
   *  @param X The point at which to evaluate the Jacobian
   *  @return A pair containing the function value and the Jacobian matrix at the given point
   *
   *  With the Sacado backend, F is evaluated once; otherwise once per component of X.
   */
  template < int n >
  std::pair< Eigen::Vector< double, n >, Eigen::Matrix< double, n, n > > dF_dX(
    const fixed_size_vector_to_vector_function_type< n >& F,
    const Eigen::Vector< double, n >&                     X )
  {
    using scalar = JacobianDual< n >;

    Eigen::Vector< scalar, n >    X_ = X.template cast< scalar >();
    Eigen::Vector< scalar, n >    F_right;
    Eigen::Vector< double, n >    F_;
    Eigen::Matrix< double, n, n > J;

    if constexpr ( isVectorDual< scalar > ) {
      for ( int j = 0; j < n; j++ )
        seedDirection( X_( j ), j );
      F_right = F( X_ );
      for ( int i = 0; i < n; i++ ) {
        F_( i ) = value( F_right( i ) );
        for ( int j = 0; j < n; j++ )
          J( i, j ) = partialDerivative( F_right( i ), j );
      }
    }
    else {
      // J_ij = d F_i / d x_j
      for ( int j = 0; j < n; j++ ) {
        seed< 1 >( X_( j ), 1.0 );
        F_right = F( X_ );
        for ( int i = 0; i < n; i++ )
          J( i, j ) = derivative< 1 >( F_right( i ) );
        seed< 1 >( X_( j ), 0.0 );
      }
      for ( int i = 0; i < n; i++ )
        F_( i ) = value( F_right( i ) );
    }

    return { F_, J };
  }

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
} // namespace Marmot::Math::Differentiation::Automatic
