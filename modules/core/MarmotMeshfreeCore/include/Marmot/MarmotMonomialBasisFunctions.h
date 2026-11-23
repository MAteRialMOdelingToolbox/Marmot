/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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

#include <Eigen/Core>
#include <vector>

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
   * @brief The complete monomial basis @f$ \boldsymbol{H}(\boldsymbol{x}) @f$ of order @f$ n @f$ in @f$ d @f$
   * dimensions, as the list of the factors of its entries.
   *
   * The entries are all monomials @f$ x_1^{\alpha_1} \cdots x_d^{\alpha_d} @f$ with
   * @f$ \alpha_1 + \dots + \alpha_d \le n @f$, ordered with the exponent of the last coordinate varying slowest
   * and that of the first coordinate fastest. For example, in 2D with @f$ n = 2 @f$:
   * @f$ \boldsymbol{H} = [1,\ x_1,\ x_1^2,\ x_2,\ x_1 x_2,\ x_2^2]^T @f$; in 3D with @f$ n = 1 @f$:
   * @f$ \boldsymbol{H} = [1,\ x_1,\ x_2,\ x_3]^T @f$. The first entry is always @f$ 1 @f$.
   *
   * The exponents @f$ \boldsymbol{\alpha} @f$ of every entry are listed once, at construction, as the factors
   * @f$ x_j^{\alpha_{kj}} @f$ of nonzero exponent of each entry; an evaluation is then the product of the powers, @f$
   * H_k = \prod_j x_j^{\alpha_{kj}} @f$, and its gradient that product with the factor of coordinate @f$ i @f$ replaced
   * by its derivative, @f$ \partial H_k / \partial x_i = \alpha_{ki}\, x_i^{\alpha_{ki} - 1} \prod_{j \neq i}
   * x_j^{\alpha_{kj}} @f$. The powers are formed by multiplication, and the factors are multiplied in the order of the
   * coordinates. Only the coordinates of nonzero exponent are factors of an entry: a factor @f$ x_j^0 = 1 @f$ would not
   * change the product, and the derivative with respect to a coordinate that is not a factor is zero. Most exponents of
   * a low-order basis are zero, e.g. all but three of the twelve of the linear basis in 3D.
   */
  class MonomialBasis {
  public:
    /**
     * @brief List the exponents of the basis.
     * @param[in] order Completeness order @f$ n \ge 0 @f$.
     * @param[in] dim   Spatial dimension @f$ d \ge 1 @f$.
     */
    MonomialBasis( int order, int dim ) : _dim( dim ), _size( computeSizeOfMonomialBasisVector( order, dim ) )
    {
      // entry k d + j: the exponent of coordinate j in monomial k
      std::vector< int > exponents( _size * dim, 0 );
      int                entry = 0;
      listExponents( order, dim, entry, exponents );

      _firstFactorOf.push_back( 0 );
      for ( int k = 0; k < _size; k++ ) {
        for ( int j = 0; j < _dim; j++ )
          if ( exponents[k * _dim + j] > 0 )
            _factors.push_back( { j, exponents[k * _dim + j] } );
        _firstFactorOf.push_back( static_cast< int >( _factors.size() ) );
      }
    }

    /// @brief The number of monomials, @f$ \binom{d + n}{d} @f$.
    int size() const { return _size; }

    /**
     * @brief Evaluate the basis.
     * @param[in]  x Coordinates @f$ \boldsymbol{x} @f$, @f$ d @f$ of them.
     * @param[out] H The basis; of size size(). Any contiguous vector, e.g. one of a fixed largest size, which needs no
     *               allocation.
     */
    void evaluate( const Eigen::Ref< const Eigen::VectorXd >& x, Eigen::Ref< Eigen::VectorXd > H ) const
    {
      for ( int k = 0; k < _size; k++ ) {
        H( k ) = 1.0;
        for ( int f = _firstFactorOf[k]; f < _firstFactorOf[k + 1]; f++ )
          H( k ) *= power( x( _factors[f].coordinate ), _factors[f].exponent );
      }
    }

    /**
     * @brief Evaluate the gradient @f$ \partial \boldsymbol{H} / \partial \boldsymbol{x} @f$ of the basis.
     * @param[in]  x     Coordinates @f$ \boldsymbol{x} @f$, @f$ d @f$ of them.
     * @param[out] dH_dx The gradient, size() @f$ \times d @f$: row @f$ k @f$ is entry @f$ k @f$ of the basis, column
     *                   @f$ i @f$ the derivative with respect to @f$ x_i @f$.
     */
    void evaluateGradient( const Eigen::Ref< const Eigen::VectorXd >& x, Eigen::Ref< Eigen::MatrixXd > dH_dx ) const
    {
      for ( int k = 0; k < _size; k++ ) {
        // zero with respect to a coordinate that is not a factor of the entry
        dH_dx.row( k ).setZero();
        // with respect to the coordinate of factor g: the product with that factor differentiated
        for ( int g = _firstFactorOf[k]; g < _firstFactorOf[k + 1]; g++ ) {
          double derivative = 1.0;
          for ( int f = _firstFactorOf[k]; f < _firstFactorOf[k + 1]; f++ ) {
            const int j     = _factors[f].coordinate;
            const int alpha = _factors[f].exponent;
            derivative *= f == g ? alpha * power( x( j ), alpha - 1 ) : power( x( j ), alpha );
          }
          dH_dx( k, _factors[g].coordinate ) = derivative;
        }
      }
    }

  private:
    /// @brief @f$ x^\alpha @f$, @f$ \alpha \ge 0 @f$, by multiplication.
    static double power( double x, int alpha )
    {
      double result = 1.0;
      for ( int p = 0; p < alpha; p++ )
        result *= x;
      return result;
    }

    /**
     * @brief List the exponents of the coordinates @f$ 1, \dots, @f$ @p dim of all monomials of total degree
     * @f$ \le @f$ @p order, from entry @p entry on: the exponent of coordinate @p dim varies slowest, the lower
     * coordinates are listed recursively.
     * @param[in] order The largest total degree.
     * @param[in] dim The number of coordinates listed.
     * @param[in,out] entry The next entry; advanced past the entries listed.
     * @param[in,out] exponents Entry @f$ k\, d + j @f$: the exponent of coordinate @f$ j @f$ in monomial @f$ k @f$.
     */
    void listExponents( int order, int dim, int& entry, std::vector< int >& exponents ) const
    {
      for ( int i = 0; i <= order; i++ ) {
        const int first = entry;
        if ( dim > 1 )
          listExponents( order - i, dim - 1, entry, exponents );
        else
          entry++;
        for ( int k = first; k < entry; k++ )
          exponents[k * _dim + dim - 1] = i;
      }
    }

    /// @brief A factor @f$ x_j^\alpha @f$, @f$ \alpha > 0 @f$, of a monomial.
    struct Factor {
      int coordinate; ///< @f$ j @f$
      int exponent;   ///< @f$ \alpha @f$
    };

    int                   _dim;     ///< the spatial dimension
    int                   _size;    ///< the number of monomials
    std::vector< Factor > _factors; ///< the factors of all monomials, monomial after monomial, in coordinate order
    /// the factors of monomial @f$ k @f$ are the entries [_firstFactorOf[k], _firstFactorOf[k + 1]) of _factors
    std::vector< int > _firstFactorOf;
  };

  /**
   * @brief Evaluate the complete monomial basis @f$ \boldsymbol{H}(\boldsymbol{x}) @f$ of order @f$ n @f$ (see
   * MonomialBasis, which lists its exponents once for many evaluations).
   *
   * @param[in]  order Completeness order @f$ n @f$.
   * @param[in]  x     Coordinates @f$ \boldsymbol{x} @f$ (the dimension @f$ d @f$ is taken from its size).
   * @param[out] res   Basis vector; must already have the size computeSizeOfMonomialBasisVector( order, d ).
   */
  inline void computeMonomialBasis( int                                        order,
                                    const Eigen::Ref< const Eigen::VectorXd >& x,
                                    Eigen::Ref< Eigen::VectorXd >              res )
  {
    MonomialBasis( order, static_cast< int >( x.size() ) ).evaluate( x, res );
  }

  /**
   * @brief Evaluate the gradient @f$ \partial \boldsymbol{H} / \partial \boldsymbol{x} @f$ of the monomial basis (see
   * MonomialBasis).
   *
   * Row @f$ k @f$ corresponds to entry @f$ k @f$ of computeMonomialBasis() (same ordering), column @f$ i @f$ to
   * @f$ \partial / \partial x_i @f$.
   *
   * @param[in]  order Completeness order @f$ n @f$.
   * @param[in]  x     Coordinates @f$ \boldsymbol{x} @f$ (the dimension @f$ d @f$ is taken from its size).
   * @param[out] res   Gradient matrix; must already have the size computeSizeOfMonomialBasisVector( order, d )
   *                   @f$ \times d @f$.
   */
  inline void computeMonomialBasisGradient( int                                        order,
                                            const Eigen::Ref< const Eigen::VectorXd >& x,
                                            Eigen::Ref< Eigen::MatrixXd >              res )
  {
    MonomialBasis( order, static_cast< int >( x.size() ) ).evaluateGradient( x, res );
  }

} // namespace Marmot::Math
