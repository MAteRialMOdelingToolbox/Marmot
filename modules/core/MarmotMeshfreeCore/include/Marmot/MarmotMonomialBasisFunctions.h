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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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

#include <Eigen/Core>
#include <cmath>

namespace Marmot::Math {

  /**
   * @brief Compute the size of the complete monomial basis of a given order.
   *
   * Computes the size of the H vector based on the completeness order and the dimension, cf. Eq. (3.67) in the
   * book by Belytschko, Chen, Hillman, i.e., the number of monomials @f$ \boldsymbol{x}^{\boldsymbol{\alpha}} @f$
   * with @f$ |\boldsymbol{\alpha}| \le n @f$ in @f$ d @f$ dimensions, @f$ \binom{d + n}{d} @f$.
   *
   * @param[in] order Completeness order @f$ n @f$ (zero for a negative order).
   * @param[in] dim   Spatial dimension @f$ d \ge 1 @f$.
   * @return The number of monomials.
   */
  inline int computeSizeOfMonomialBasisVector( int order, int dim )
  {
    // compute the size of the H vector based on the completeness order and the dimension
    // for Eq. (3.67) in the book by Belytschko, Chen, Hillman.

    int size = 0;
    for ( int i = 0; i <= order; i++ ) {
      if ( dim == 1 )
        size += 1;
      else
        size += computeSizeOfMonomialBasisVector( order - i, dim - 1 );
    }
    return size;
  }

  /**
   * @brief Recursive kernel of computeMonomialBasis().
   *
   * Multiplies the entries starting at @p idxEnd by the powers @f$ x_k^i @f$ of the coordinates
   * @f$ k \le @f$ @p dim, for all exponent combinations of total degree @f$ \le @f$ @p order. The loop over the
   * exponent of coordinate @p dim is the outermost one, the recursion handles the lower coordinates.
   *
   * @param[in]     order  Remaining total degree.
   * @param[in]     x      Coordinates.
   * @param[in,out] res    Basis vector, initialized with ones by the caller.
   * @param[in]     idxEnd First entry of @p res treated by this call.
   * @param[in]     dim    Number of (leading) coordinates treated by this call.
   * @return One past the last entry of @p res treated by this call.
   */
  inline int _computeMonomialBasisRecursion( int                    order,
                                             const Eigen::VectorXd& x,
                                             Eigen::VectorXd&       res,
                                             int                    idxEnd,
                                             int                    dim )
  {
    for ( int i = 0; i <= order; i++ ) {

      const int idxStart = idxEnd;
      if ( dim > 1 )
        idxEnd = _computeMonomialBasisRecursion( order - i, x, res, idxEnd, dim - 1 );
      else {
        idxEnd++;
      }

      for ( int idx = idxStart; idx < idxEnd; idx++ )
        res( idx ) *= std::pow( x[dim - 1], i );
    }
    return idxEnd;
  }

  /**
   * @brief Recursive kernel of computeMonomialBasisGradient().
   *
   * Same traversal as _computeMonomialBasisRecursion(). By the product rule, column @f$ k @f$ of an entry is
   * multiplied by @f$ \partial x_k^i / \partial x_k = i\, x_k^{i-1} @f$ for its own coordinate and by
   * @f$ x_j^i @f$ for all other coordinates @f$ j \neq k @f$.
   *
   * @param[in]     order  Remaining total degree.
   * @param[in]     x      Coordinates.
   * @param[in,out] res    Gradient matrix (basis size @f$ \times d @f$), initialized with ones by the caller.
   * @param[in]     idxEnd First row of @p res treated by this call.
   * @param[in]     dim    Number of (leading) coordinates treated by this call.
   * @return One past the last row of @p res treated by this call.
   */
  inline int _computeMonomialBasisGradientRecursion( int                    order,
                                                     const Eigen::VectorXd& x,
                                                     Eigen::MatrixXd&       res,
                                                     int                    idxEnd,
                                                     int                    dim )
  {
    // the basis has the same layout as in _computeMonomialBasisRecursion; by the product rule, column k of the
    // gradient is the product of the factors x_d^i of all dimensions d, with x_k^i replaced by its derivative
    for ( int i = 0; i <= order; i++ ) {

      const int idxStart = idxEnd;

      if ( dim > 1 )
        idxEnd = _computeMonomialBasisGradientRecursion( order - i, x, res, idxEnd, dim - 1 );
      else {
        idxEnd++;
      }

      const double factor           = std::pow( x[dim - 1], i );
      const double factorDerivative = i > 0 ? i * std::pow( x[dim - 1], i - 1 ) : 0.0;

      for ( int idx = idxStart; idx < idxEnd; idx++ )
        for ( int k = 0; k < res.cols(); k++ )
          res( idx, k ) *= k == dim - 1 ? factorDerivative : factor;
    }
    return idxEnd;
  }

  /**
   * @brief Evaluate the complete monomial basis @f$ \boldsymbol{H}(\boldsymbol{x}) @f$ of order @f$ n @f$.
   *
   * The entries are all monomials @f$ x_1^{\alpha_1} \cdots x_d^{\alpha_d} @f$ with
   * @f$ \alpha_1 + \dots + \alpha_d \le n @f$, ordered with the exponent of the last coordinate varying slowest
   * and that of the first coordinate fastest. For example, in 2D with @f$ n = 2 @f$:
   * @f$ \boldsymbol{H} = [1,\ x_1,\ x_1^2,\ x_2,\ x_1 x_2,\ x_2^2]^T @f$; in 3D with @f$ n = 1 @f$:
   * @f$ \boldsymbol{H} = [1,\ x_1,\ x_2,\ x_3]^T @f$. The first entry is always @f$ 1 @f$.
   *
   * @param[in]  order Completeness order @f$ n @f$.
   * @param[in]  x     Coordinates @f$ \boldsymbol{x} @f$ (the dimension @f$ d @f$ is taken from its size).
   * @param[out] res   Basis vector; must already have the size computeSizeOfMonomialBasisVector( order, d ).
   */
  inline void computeMonomialBasis( int order, const Eigen::VectorXd& x, Eigen::VectorXd& res )
  {
    res.setOnes();
    _computeMonomialBasisRecursion( order, x, res, 0, x.size() );
  }

  /**
   * @brief Evaluate the gradient @f$ \partial \boldsymbol{H} / \partial \boldsymbol{x} @f$ of the monomial basis.
   *
   * Row @f$ k @f$ corresponds to entry @f$ k @f$ of computeMonomialBasis() (same ordering), column @f$ i @f$ to
   * @f$ \partial / \partial x_i @f$.
   *
   * @param[in]  order Completeness order @f$ n @f$.
   * @param[in]  x     Coordinates @f$ \boldsymbol{x} @f$ (the dimension @f$ d @f$ is taken from its size).
   * @param[out] res   Gradient matrix; must already have the size computeSizeOfMonomialBasisVector( order, d )
   *                   @f$ \times d @f$.
   */
  inline void computeMonomialBasisGradient( int order, const Eigen::VectorXd& x, Eigen::MatrixXd& res )
  {
    res.setOnes();
    _computeMonomialBasisGradientRecursion( order, x, res, 0, x.size() );
  }

} // namespace Marmot::Math
