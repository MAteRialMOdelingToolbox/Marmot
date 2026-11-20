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
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief The rows of a particle's tangent that belong to one of its nodes, factored by the chain rule through the
   * increment the particle interpolates.
   *
   * @details A particle sees the increments @f$ \Delta\boldsymbol{Q}_B @f$ of its nodes only through the increment it
   * interpolates at its evaluation point: its value, its gradient and, for a stabilized particle, its second gradient,
   * @f[
   *   \Delta\boldsymbol{q} = \sum_B N_B\, \Delta\boldsymbol{Q}_B, \qquad
   *   \frac{\partial \Delta\boldsymbol{q}}{\partial Y_l} = \sum_B \frac{\partial N_B}{\partial Y_l}\,
   *   \Delta\boldsymbol{Q}_B, \qquad
   *   \frac{\partial^2 \Delta\boldsymbol{q}}{\partial Y_l \partial Y_m} = \sum_B
   *   \frac{\partial^2 N_B}{\partial Y_l \partial Y_m}\, \Delta\boldsymbol{Q}_B .
   * @f]
   * By the chain rule, the block of the tangent that couples the force @f$ \boldsymbol{f}_A @f$ of node @f$ A @f$
   * with the increment of node @f$ B @f$ is
   * @f[
   *   \boldsymbol{K}_{AB} = \frac{\partial \boldsymbol{f}_A}{\partial \Delta\boldsymbol{q}}\, N_B
   *   + \sum_l \frac{\partial \boldsymbol{f}_A}{\partial (\partial \Delta\boldsymbol{q} / \partial Y_l)}\,
   *   \frac{\partial N_B}{\partial Y_l}
   *   + \sum_{l,m} \frac{\partial \boldsymbol{f}_A}{\partial (\partial^2 \Delta\boldsymbol{q} / \partial Y_l
   *   \partial Y_m)}\, \frac{\partial^2 N_B}{\partial Y_l \partial Y_m} .
   * @f]
   * The three derivatives depend on node @f$ A @f$ (through its test function) and on the state of the particle, but
   * not on @f$ B @f$: a particle computes them once per node, and every block of the node's rows follows from the
   * shape functions of the other node (addTangentRowsOfNode()). The same factors give a cheap bound of the absolute
   * value of the rows (AbsoluteTangentBound).
   *
   * Each derivative is a matrix of nDofPerNode x nDofPerNode: row @f$ i @f$ is the component of the force of the
   * node, column @f$ k @f$ the component of the interpolated increment (the fields of the node, in its dof order).
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   */
  template < int nDim, int nDofPerNode >
  struct TangentRowsOfNode {
    /// @brief A derivative of the force of the node with respect to a component of the interpolated increment.
    using Block = Eigen::Matrix< double, nDofPerNode, nDofPerNode >;

    /// @f$ \partial \boldsymbol{f}_A / \partial \Delta\boldsymbol{q} @f$.
    Block df_dIncrement;
    /// Entry @f$ l @f$: @f$ \partial \boldsymbol{f}_A / \partial (\partial \Delta\boldsymbol{q} / \partial Y_l) @f$.
    std::array< Block, nDim > df_dIncrementGradient;
    /// Entry @f$ l\, n_{\text{dim}} + m @f$: @f$ \partial \boldsymbol{f}_A / \partial (\partial^2
    /// \Delta\boldsymbol{q} / \partial Y_l \partial Y_m) @f$; zero for a particle without second gradients.
    std::array< Block, nDim * nDim > df_dIncrementSecondGradient;

    /// @brief Set all derivatives to zero.
    void setZero()
    {
      df_dIncrement.setZero();
      for ( auto& block : df_dIncrementGradient )
        block.setZero();
      for ( auto& block : df_dIncrementSecondGradient )
        block.setZero();
    }
  };

  /**
   * @brief Add the rows of node @p A to the tangent of a particle: @f$ \boldsymbol{K}_{AB} @f$ for all nodes
   * @f$ B @f$, from the factors of TangentRowsOfNode and the trial shape functions at the evaluation point.
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   * @param[in] rows The derivatives of the force of node @p A with respect to the interpolated increment.
   * @param[in] A The node.
   * @param[in] N The trial shape functions @f$ N_B @f$, 1 x nNodes.
   * @param[in] dN_dY Their gradients @f$ \partial N_B / \partial Y_l @f$, nDim x nNodes.
   * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes, entry @f$ (l\, n_{\text{dim}} + m, B) @f$ is
   * @f$ \partial^2 N_B / \partial Y_l \partial Y_m @f$; nullptr for a particle without second gradients.
   * @param[in,out] K The tangent, column-major, (nNodes nDofPerNode) x (nNodes nDofPerNode); the rows are added.
   */
  template < int nDim, int nDofPerNode >
  void addTangentRowsOfNode( const TangentRowsOfNode< nDim, nDofPerNode >& rows,
                             int                                           A,
                             const Eigen::MatrixXd&                        N,
                             const Eigen::MatrixXd&                        dN_dY,
                             const Eigen::MatrixXd*                        d2N_dYdY,
                             Eigen::Ref< Eigen::MatrixXd >                 K )
  {
    const int idxA = nDofPerNode * A;
    for ( int B = 0; B < N.size(); B++ ) {
      typename TangentRowsOfNode< nDim, nDofPerNode >::Block K_AB = rows.df_dIncrement * N( B );
      for ( int l = 0; l < nDim; l++ )
        K_AB += rows.df_dIncrementGradient[l] * dN_dY( l, B );
      if ( d2N_dYdY )
        for ( int lm = 0; lm < nDim * nDim; lm++ )
          K_AB += rows.df_dIncrementSecondGradient[lm] * ( *d2N_dYdY )( lm, B );
      K.template block< nDofPerNode, nDofPerNode >( idxA, nDofPerNode * B ) += K_AB;
    }
  }

  /**
   * @brief An upper bound of the absolute value of a particle's tangent times nonnegative weights,
   * @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$, from the factors of TangentRowsOfNode, without the tangent.
   *
   * @details What a stable time increment needs from a tangent is a bound of its spectral radius, e.g. the
   * Collatz-Wielandt bound, @f$ \max_i (|\boldsymbol{K}|\boldsymbol{w})_i / w_i @f$ for positive weights
   * @f$ \boldsymbol{w} @f$. The tangent of a particle with @f$ n @f$ nodes has @f$ n^2 @f$ blocks; this bound of
   * @f$ |\boldsymbol{K}|\boldsymbol{W} @f$ costs @f$ O(n) @f$: the rows of node @f$ A @f$ are linear in the shape
   * functions of the other node @f$ B @f$ (see TangentRowsOfNode), and the sum over @f$ B @f$ is taken from a few
   * weighted moments of the shape functions, computed once per particle. Row @f$ i @f$ of node @f$ A @f$, column
   * component @f$ k @f$, with the weights @f$ w_B @f$ of component @f$ k @f$ of the nodes and the factors
   * @f$ a = \partial f_{Ai}/\partial \Delta q_k @f$, @f$ \boldsymbol{b} @f$ (its gradient factors) and
   * @f$ \boldsymbol{h} @f$ (its second-gradient factors):
   * @f[
   *   \sum_B |K_{Ai,Bk}|\, w_B \le |a| \sum_B w_B |N_B|
   *   + \sqrt{ \sum_B w_B \|\boldsymbol{g}_B\| \; \boldsymbol{b}^T \Big( \sum_B
   *     \frac{w_B}{\|\boldsymbol{g}_B\|} \boldsymbol{g}_B \boldsymbol{g}_B^T \Big) \boldsymbol{b} }
   *   + (\text{the same for } \boldsymbol{h} \text{ and the second gradients}),
   * @f]
   * with @f$ \boldsymbol{g}_B = \partial N_B / \partial \boldsymbol{Y} @f$. The first term is exact; the others are
   * the Cauchy-Schwarz inequality over the nodes, @f$ \sum_B w_B |\boldsymbol{b}\cdot\boldsymbol{g}_B| \le
   * \sqrt{\sum_B w_B \|\boldsymbol{g}_B\|} \sqrt{\sum_B w_B (\boldsymbol{b}\cdot\boldsymbol{g}_B)^2 /
   * \|\boldsymbol{g}_B\|} @f$, which keeps the directions of the gradients: it is exact if all gradients are
   * parallel, and @f$ 2/\sqrt{3} @f$ times the exact sum (@f$ \pi/(2\sqrt{2}) @f$ in 2D) if their directions are
   * spread evenly. The terms are added by the triangle inequality. The bound holds for every tangent of the form of
   * TangentRowsOfNode, symmetric or not; the second gradients of the shape functions must be symmetric (their
   * symmetric components are used).
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   */
  template < int nDim, int nDofPerNode >
  class AbsoluteTangentBound {
  public:
    /// The number of independent components of a symmetric second gradient.
    static constexpr int nSym = nDim * ( nDim + 1 ) / 2;

    /**
     * @brief Compute the weighted moments of the shape functions.
     * @param[in] columnWeights The nonnegative weights @f$ \boldsymbol{W} @f$, column-major, (nNodes nDofPerNode) x
     * nColumns, in the dof order of the particle.
     * @param[in] nColumns The number of columns of @p columnWeights.
     * @param[in] N The trial shape functions @f$ N_B @f$, 1 x nNodes.
     * @param[in] dN_dY Their gradients, nDim x nNodes.
     * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes (see addTangentRowsOfNode()); nullptr for a
     * particle without second gradients.
     * @param[in] nComponentsWithSecondGradients The components of a node whose columns have second-gradient factors
     * are its first nComponentsWithSecondGradients (e.g. the displacement, on which a stabilization acts, but not a
     * nonlocal field): the moments of the second gradients are formed for these columns only. addRowsOfNode() refuses
     * a second-gradient factor in another column.
     */
    AbsoluteTangentBound( const double*          columnWeights,
                          int                    nColumns,
                          const Eigen::MatrixXd& N,
                          const Eigen::MatrixXd& dN_dY,
                          const Eigen::MatrixXd* d2N_dYdY,
                          int                    nComponentsWithSecondGradients = nDofPerNode )
      : _nDof( nDofPerNode * static_cast< int >( N.size() ) ),
        _nColumns( nColumns ),
        _withSecondGradients( d2N_dYdY != nullptr ),
        _nComponentsWithSecondGradients( nComponentsWithSecondGradients ),
        _moments( nColumns * nDofPerNode )
    {
      const int nNodes = static_cast< int >( N.size() );
      for ( int c = 0; c < nColumns; c++ )
        for ( int k = 0; k < nDofPerNode; k++ ) {
          Moments& m = _moments[c * nDofPerNode + k];
          for ( int B = 0; B < nNodes; B++ ) {
            const double w_B = columnWeights[c * _nDof + nDofPerNode * B + k];
            if ( w_B == 0.0 )
              continue;
            m.active = true;
            m.sumOfValues += w_B * std::abs( N( B ) );
            const Eigen::Matrix< double, nDim, 1 > g_B = dN_dY.col( B );
            addMoments( g_B, w_B, m.sumOfGradientNorms, m.gradientMoments );
            if ( _withSecondGradients && k < _nComponentsWithSecondGradients )
              addMoments( symmetricComponents( d2N_dYdY->col( B ) ),
                          w_B,
                          m.sumOfSecondGradientNorms,
                          m.secondGradientMoments );
          }
        }
    }

    /**
     * @brief Add the bound of the rows of node @p A: @f$ (|\boldsymbol{K}|\boldsymbol{W})_{Ai,c} @f$ for all
     * components @f$ i @f$ and columns @f$ c @f$.
     * @param[in] rows The derivatives of the force of node @p A with respect to the interpolated increment.
     * @param[in] A The node.
     * @param[in,out] bound The bound, column-major, (nNodes nDofPerNode) x nColumns; the bound of the rows is added.
     */
    void addRowsOfNode( const TangentRowsOfNode< nDim, nDofPerNode >& rows, int A, double* bound ) const
    {
      for ( int i = 0; i < nDofPerNode; i++ )
        for ( int k = 0; k < nDofPerNode; k++ ) {
          const double                     a = rows.df_dIncrement( i, k );
          Eigen::Matrix< double, nDim, 1 > b;
          for ( int l = 0; l < nDim; l++ )
            b( l ) = rows.df_dIncrementGradient[l]( i, k );
          Eigen::Matrix< double, nSym, 1 > h = Eigen::Matrix< double, nSym, 1 >::Zero();
          if ( _withSecondGradients )
            for ( int l = 0, s = 0; l < nDim; l++ )
              for ( int m = l; m < nDim; m++, s++ )
                h( s ) = rows.df_dIncrementSecondGradient[l * nDim + m]( i, k ) +
                         ( m != l ? rows.df_dIncrementSecondGradient[m * nDim + l]( i, k ) : 0.0 );

          // Most factors of a node's rows are structurally zero (e.g. a displacement row has no gradient factor for a
          // nonlocal column, and only the stabilization has second-gradient factors): a zero factor adds nothing.
          const bool withValueTerm          = a != 0.0;
          const bool withGradientTerm       = ( b.array() != 0.0 ).any();
          const bool withSecondGradientTerm = ( h.array() != 0.0 ).any();
          if ( !withValueTerm && !withGradientTerm && !withSecondGradientTerm )
            continue;
          if ( withSecondGradientTerm && k >= _nComponentsWithSecondGradients )
            throw std::logic_error( "AbsoluteTangentBound: a second-gradient factor in the column of a component "
                                    "without second gradients, whose moments are not formed" );

          for ( int c = 0; c < _nColumns; c++ ) {
            const Moments& m = _moments[c * nDofPerNode + k];
            if ( !m.active )
              continue;
            double sum = std::abs( a ) * m.sumOfValues;
            if ( withGradientTerm )
              sum += std::sqrt( std::max( 0.0, m.sumOfGradientNorms * b.dot( m.gradientMoments * b ) ) );
            if ( withSecondGradientTerm )
              sum += std::sqrt( std::max( 0.0, m.sumOfSecondGradientNorms * h.dot( m.secondGradientMoments * h ) ) );
            bound[c * _nDof + nDofPerNode * A + i] += sum;
          }
        }
    }

  private:
    /// @brief The weighted moments of the shape functions for one column and one component of the nodes.
    struct Moments {
      bool   active             = false; ///< some weight of the column and component is positive
      double sumOfValues        = 0.0;   ///< @f$ \sum_B w_B |N_B| @f$
      double sumOfGradientNorms = 0.0;   ///< @f$ \sum_B w_B \|\boldsymbol{g}_B\| @f$
      /// @f$ \sum_B w_B \boldsymbol{g}_B \boldsymbol{g}_B^T / \|\boldsymbol{g}_B\| @f$
      Eigen::Matrix< double, nDim, nDim > gradientMoments          = Eigen::Matrix< double, nDim, nDim >::Zero();
      double                              sumOfSecondGradientNorms = 0.0; ///< the same for the second gradients
      /// the same for the symmetric components of the second gradients
      Eigen::Matrix< double, nSym, nSym > secondGradientMoments = Eigen::Matrix< double, nSym, nSym >::Zero();
    };

    /**
     * @brief Add the moments of one node: @f$ w \|\boldsymbol{g}\| @f$ and @f$ w\, \boldsymbol{g}
     * \boldsymbol{g}^T / \|\boldsymbol{g}\| @f$ (nothing for @f$ \boldsymbol{g} = \boldsymbol{0} @f$).
     * @tparam n The size of @p g.
     * @param[in] g The gradient (or the symmetric components of the second gradient) of the node.
     * @param[in] w The weight of the node.
     * @param[in,out] sumOfNorms The sum of the weighted norms.
     * @param[in,out] moments The sum of the weighted outer products.
     */
    template < int n >
    static void addMoments( const Eigen::Matrix< double, n, 1 >& g,
                            double                               w,
                            double&                              sumOfNorms,
                            Eigen::Matrix< double, n, n >&       moments )
    {
      const double norm = g.norm();
      if ( norm == 0.0 )
        return;
      sumOfNorms += w * norm;
      moments += ( w / norm ) * g * g.transpose();
    }

    /**
     * @brief The independent components @f$ (l, m), l \le m @f$ of a symmetric second gradient.
     * @param[in] d2N_dYdY The second gradient of one node, entry @f$ l\, n_{\text{dim}} + m @f$.
     * @return The nSym components, in the order @f$ (0,0), (0,1), \ldots, (1,1), \ldots @f$.
     */
    static Eigen::Matrix< double, nSym, 1 > symmetricComponents( const Eigen::Ref< const Eigen::VectorXd >& d2N_dYdY )
    {
      Eigen::Matrix< double, nSym, 1 > components;
      for ( int l = 0, s = 0; l < nDim; l++ )
        for ( int m = l; m < nDim; m++, s++ )
          components( s ) = d2N_dYdY( l * nDim + m );
      return components;
    }

    const int              _nDof;                           ///< the number of dofs of the particle
    const int              _nColumns;                       ///< the number of columns of the weights
    const bool             _withSecondGradients;            ///< the particle has second gradients
    const int              _nComponentsWithSecondGradients; ///< the leading components with second-gradient factors
    std::vector< Moments > _moments;                        ///< entry @f$ c\, n_{\text{dof/node}} + k @f$
  };

} // namespace Marmot::Discretizations::Meshfree
