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
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include <functional>

namespace Marmot {
  namespace NumericalAlgorithms::Differentiation {

    /** @typedef scalar_to_scalar_function_type
     *  @brief Type definition for a function that maps a double to another double.
     */
    using scalar_to_scalar_function_type = std::function< double( const double x ) >;

    /** @typedef vector_to_vector_function_type
     *  @brief Type definition for a function that maps an Eigen vector to another Eigen vector.
     */
    using vector_to_vector_function_type = std::function< Eigen::VectorXd( const Eigen::VectorXd& X ) >;

    /** @brief Approximates the first derivative of a scalar function f at point x using the forward difference method.
     *  @param f The scalar function for which the derivative is to be approximated.
     *  @param x The point at which the derivative is to be approximated.
     *  @return The approximated first derivative of f at point x.
     *
     *  It actually computes
     *  \f[
     *   f'(x) \approx \frac{f(x + h) - f(x)}{h}
     *   \f]
     *   with \f$ h = \sqrt{\epsilon}\, \max(1,|x|) \f$ and \f$ \epsilon \f$ being the machine precision.
     */
    double forwardDifference( const scalar_to_scalar_function_type& f, const double x );

    /** @brief Approximates the first derivative of a scalar function f at point x using the central difference method.
     *  @param f The scalar function for which the derivative is to be approximated.
     *  @param x The point at which the derivative is to be approximated.
     *  @return The approximated first derivative of f at point x.
     *
     *  It actually approximates the first derivative by
     *  \f[
     *   f'(x) \approx \frac{f(x + h) - f(x - h)}{2h}
     *   \f]
     *   with \f$ h = \sqrt[3]{\epsilon}\,\max( 1, |x|)\f$ and \f$ \epsilon \f$ being the machine precision.
     */
    double centralDifference( const scalar_to_scalar_function_type& f, const double x );

    /** @brief Approximates the Jacobian matrix of a vector function F at point X using the forward difference method.
     *  @param F The vector function for which the Jacobian matrix is to be approximated.
     *  @param X The point at which the Jacobian matrix is to be approximated.
     *  @return The approximated Jacobian matrix of F at point X.
     *
     *  It actually computes
     *  \f[
     *   J_{ij} \approx \frac{F_i(X + h e_j) - F_i(X)}{h}
     *   \f]
     *   with \f$ h = \sqrt{\epsilon} (1 + |X_j|) \f$, \f$ e_j \f$ being the unit vector in the j-th direction, and \f$
     * \epsilon \f$ being the machine precision.
     */
    Eigen::MatrixXd forwardDifference( const vector_to_vector_function_type& F, const Eigen::VectorXd& X );

    /** @brief Approximates the Jacobian matrix of a vector function F at point X using the central difference method.
     *  @param F The vector function for which the Jacobian matrix is to be approximated.
     *  @param X The point at which the Jacobian matrix is to be approximated.
     *  @return The approximated Jacobian matrix of F at point X.
     *
     *  It actually approximates the Jocobian matrix by
     *  \f[
     *   J_{ij} \approx \frac{F_i(\boldsymbol{X} + h e_j) - F_i(\boldsymbol{X} - h e_j)}{2h}
     *   \f]
     *   with \f$ h = \sqrt[3]{\epsilon} \max(1, ||\boldsymbol{X}||) \f$, \f$ e_j \f$ being the unit vector in the j-th
     * direction, and \f$ \epsilon \f$ being the machine precision.
     */
    Eigen::MatrixXd centralDifference( const vector_to_vector_function_type& F, const Eigen::VectorXd& X );

    namespace Complex {

      /// A variable representing \f$ 0 + 1i \f$
      const static std::complex< double > imaginaryUnit = { 0, 1 };
      /// A variable representing \f$ 1 + 1i \f$
      const static std::complex< double > complexUnit = { 1, 1 };
      /// A variable representing \f$ \frac{\sqrt{2}}{2} (1+1i) \f$
      const static std::complex< double > i_ = Marmot::Constants::sqrt2 / 2. * complexUnit;

      /** @typedef scalar_to_scalar_function_type
       *  @brief Type definition for a function that maps a complex double to another complex double.
       */
      using scalar_to_scalar_function_type = std::function< complexDouble( const complexDouble x ) >;

      /** @typedef vector_to_vector_function_type
       *  @brief Type definition for a function that maps an Eigen complex vector to another Eigen complex vector.
       */
      using vector_to_vector_function_type = std::function< Eigen::VectorXcd( const Eigen::VectorXcd& X ) >;

      /** @brief Approximates the first derivative of a scalar function f at point x using the complex step method.
       *  @param f The scalar function for which the derivative is to be approximated.
       *  @param x The point at which the derivative is to be approximated.
       *  @return The approximated first derivative of f at point x.
       *
       *  It actually computes
       *  \f[
       *   f'(x) \approx \frac{\text{Im}(f(x + ih))}{h}
       *   \f]
       *   with \f$ h = 10^{-20} \f$.
       *
       *   @note This method is highly accurate and does not suffer from subtractive cancellation errors,
       *   making it suitable for functions where high precision is required.
       *
       */
      double forwardDifference( const scalar_to_scalar_function_type& f, const double x );

      /** @brief Approximates the Jacobian matrix of a vector function F at point X using the complex step method.
       *  @param F The vector function for which the Jacobian matrix is to be approximated.
       *  @param X The point at which the Jacobian matrix is to be approximated.
       *  @return A tuple containing the function value at X and the approximated Jacobian matrix of F at point X.
       *
       *  It actually computes
       *  \f[
       *   J_{ij} \approx \frac{\text{Im}(F_i(X + ih e_j))}{h}
       *   \f]
       *   with \f$ h = 10^{-20} \f$, and \f$ e_j \f$ being the unit vector in the j-th direction.
       *
       *   @note This method is highly accurate and does not suffer from subtractive cancellation errors,
       *   making it suitable for functions where high precision is required.
       *
       */
      std::tuple< Eigen::VectorXd, Eigen::MatrixXd > forwardDifference( const vector_to_vector_function_type& F,
                                                                        const Eigen::VectorXd&                X );

      /** @brief Approximates the Jacobian matrix of a vector function F at point X using the complex step method with
       * central differences.
       *  @param F The vector function for which the Jacobian matrix is to be approximated.
       *  @param X The point at which the Jacobian matrix is to be approximated.
       *  @return The approximated Jacobian matrix of F at point X.
       *
       *  It actually approximates the Jacobian matrix by
       *  \f[
       *  J_{ij} \approx \frac{\text{Im}(F_i(\boldsymbol{X} + Ih e_j) - F_i(\boldsymbol{X} - Ih e_j))}{2h}
       *  \f]
       *  with  \f$ I = \sqrt{2}/2( 1+1i )\f$,
       *  \f$ h = \sqrt[3]{\epsilon} \max(1, ||\boldsymbol{X}||) \f$,
       *  and \f$ e_j \f$ being the unit vector in the j-th direction.
       *
       *  @note The formula can be found in Lai et al. (2005) Equ. 19.
       */
      Eigen::MatrixXd centralDifference( const vector_to_vector_function_type& F, const Eigen::VectorXd& X );

      /** @brief Approximates the Jacobian matrix of a vector function F at point X using a fourth-order accurate
       * complex step method.
       * @param F The vector function for which the Jacobian matrix is to be approximated.
       * @param X The point at which the Jacobian matrix is to be approximated.
       * @return The approximated Jacobian matrix of F at point X.
       *
       * It actually approximates the Jacobian matrix by
       * \f[
       * J_{ij} \approx \frac{\text{Im}(8 [
       * F_i(\boldsymbol{X} + I/2h e_j) -
       * F_i(\boldsymbol{X} - I/2h e_j)] -
       * [F_i(\boldsymbol{X} + Ih e_j) +
       * F_i(\boldsymbol{X} - Ih e_j)])}{3\sqrt{2}h}
       * \f]
       * with  \f$ I = \sqrt{2}/2( 1+1i )\f$,
       * \f$ h = \sqrt{\epsilon} \max(1, ||\boldsymbol{X}||) \f$,
       * and \f$ e_j \f$ being the unit vector in the j-th direction.
       *
       * @note The formula can be found in Lai et al. (2005) Equ. 24.
       */
      Eigen::MatrixXd fourthOrderAccurateDerivative( const vector_to_vector_function_type& F,
                                                     const Eigen::VectorXd&                X );

    } // namespace Complex
  }   // namespace NumericalAlgorithms::Differentiation
} // namespace Marmot
