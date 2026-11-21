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

#include "Marmot/MarmotFastorTensorBasics.h"
#include <Eigen/Core>
#include <Fastor/Fastor.h>
#include <algorithm>
#include <array>
#include <cmath>
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
   * @f$ \partial^2 N_B / \partial Y_l \partial Y_m @f$; empty for a particle without second gradients.
   * @param[in,out] K The tangent, column-major, (nNodes nDofPerNode) x (nNodes nDofPerNode); the rows are added.
   */
  template < int nDim, int nDofPerNode >
  void addTangentRowsOfNode( const TangentRowsOfNode< nDim, nDofPerNode >& rows,
                             int                                           A,
                             const Eigen::MatrixXd&                        N,
                             const Eigen::MatrixXd&                        dN_dY,
                             const Eigen::MatrixXd&                        d2N_dYdY,
                             Eigen::Ref< Eigen::MatrixXd >                 K )
  {
    const int idxA = nDofPerNode * A;
    for ( int B = 0; B < N.size(); B++ ) {
      typename TangentRowsOfNode< nDim, nDofPerNode >::Block K_AB = rows.df_dIncrement * N( B );
      for ( int l = 0; l < nDim; l++ )
        K_AB += rows.df_dIncrementGradient[l] * dN_dY( l, B );
      if ( d2N_dYdY.size() )
        for ( int lm = 0; lm < nDim * nDim; lm++ )
          K_AB += rows.df_dIncrementSecondGradient[lm] * d2N_dYdY( lm, B );
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
     * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes (see addTangentRowsOfNode()); empty for a particle
     * without second gradients.
     */
    AbsoluteTangentBound( const double*          columnWeights,
                          int                    nColumns,
                          const Eigen::MatrixXd& N,
                          const Eigen::MatrixXd& dN_dY,
                          const Eigen::MatrixXd& d2N_dYdY )
      : _nDof( nDofPerNode * static_cast< int >( N.size() ) ),
        _nColumns( nColumns ),
        _withSecondGradients( d2N_dYdY.size() != 0 ),
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
            if ( _withSecondGradients )
              addMoments( symmetricComponents( d2N_dYdY.col( B ) ),
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

    const int              _nDof;                ///< the number of dofs of the particle
    const int              _nColumns;            ///< the number of columns of the weights
    const bool             _withSecondGradients; ///< the particle has second gradients
    std::vector< Moments > _moments;             ///< entry @f$ c\, n_{\text{dof/node}} + k @f$
  };

  /**
   * @brief Add the bound of @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$ of a particle, node by node
   * (AbsoluteTangentBound): what Particle::computeBoundOfAbsoluteStiffnessTimes() gives, once the particle's material
   * is evaluated at the state the bound is for.
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   * @tparam RowsOfNode  A callable `(int A, TangentRowsOfNode& rows)` that adds the derivatives of node @p A to
   * the (zeroed) @p rows, e.g. a particle's addStiffnessRowsOfNode().
   * @param[in] columnWeights The nonnegative weights, see AbsoluteTangentBound.
   * @param[in] nColumns The number of columns of @p columnWeights.
   * @param[in] N The trial shape functions, 1 x nNodes.
   * @param[in] dN_dY Their gradients, nDim x nNodes.
   * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes; empty for a particle without second gradients.
   * @param[in] addRowsOfNode The derivatives of the force of a node.
   * @param[in,out] bound The bound, column-major, (nNodes nDofPerNode) x nColumns; the bound of the rows is added.
   */
  template < int nDim, int nDofPerNode, typename RowsOfNode >
  void addBoundOfAbsoluteTangentTimesWeights( const double*          columnWeights,
                                              int                    nColumns,
                                              const Eigen::MatrixXd& N,
                                              const Eigen::MatrixXd& dN_dY,
                                              const Eigen::MatrixXd& d2N_dYdY,
                                              const RowsOfNode&      addRowsOfNode,
                                              double*                bound )
  {
    const AbsoluteTangentBound< nDim, nDofPerNode > absoluteTangent( columnWeights, nColumns, N, dN_dY, d2N_dYdY );
    TangentRowsOfNode< nDim, nDofPerNode >          rows;
    for ( int A = 0; A < N.size(); A++ ) {
      rows.setZero();
      addRowsOfNode( A, rows );
      absoluteTangent.addRowsOfNode( rows, A, bound );
    }
  }

  /**
   * @brief The NSNI stabilization of the displacement rows of a particle, shared by the displacement and the
   * gradient-enhanced particles (whose dofs of a node start with the displacement).
   *
   * @details The stabilization adds to the force of node @f$ A @f$ the term @f$ D_{iK}\, \partial \tau_{ij} /
   * \partial Y_K @f$, with the moment-weighted second gradient of the test function
   * @f[
   *   D_{iK} = \frac{\partial Y_j}{\partial x_i}\, \frac{\partial^2 T_A}{\partial Y_j \partial Y_L}\,
   *   \frac{M^Y_{LK}}{J_Y} ,
   * @f]
   * the second moments @f$ \boldsymbol{M}^Y @f$ of the particle domain and @f$ J_Y = \det \partial\boldsymbol{Y} /
   * \partial\boldsymbol{X} @f$. The second gradients of the shape functions are the smoothed ones of the particle
   * domain, @f$ \partial^2 N_B / \partial Y_I \partial Y_J @f$ from the face sums `d2N_dYdY[J](I, B)` of
   * @f$ \partial N_B / \partial Y_I\, n_J\, dA / V @f$, symmetrized.
   */
  namespace NSNI {

    /**
     * @brief The symmetrized smoothed second gradient of the shape function of node @p B.
     * @tparam nDim Spatial dimension.
     * @param[in] d2N_dYdY The face sums: `d2N_dYdY[J](I, B)` is the sum of @f$ \partial N_B / \partial Y_I\, n_J\,
     * dA / V @f$.
     * @param[in] B The node.
     * @return @f$ \tfrac{1}{2} (\text{d2N\_dYdY}[J](I, B) + \text{d2N\_dYdY}[I](J, B)) @f$ at @f$ (I, J) @f$.
     */
    template < int nDim >
    Fastor::Tensor< double, nDim, nDim > secondGradientOfNode( const std::array< Eigen::MatrixXd, nDim >& d2N_dYdY,
                                                               int                                        B )
    {
      Fastor::Tensor< double, nDim, nDim > d2N_B_dYdY( 0.0 );
      for ( int J = 0; J < nDim; J++ )
        for ( int I = 0; I < nDim; I++ )
          d2N_B_dYdY( I, J ) += 0.5 * ( d2N_dYdY[J]( I, B ) + d2N_dYdY[I]( J, B ) );
      return d2N_B_dYdY;
    }

    /**
     * @brief The symmetrized smoothed second gradients of all nodes, in the layout of addTangentRowsOfNode().
     * @tparam nDim Spatial dimension.
     * @param[in] d2N_dYdY The face sums, see secondGradientOfNode().
     * @param[in] nNodes The number of nodes.
     * @return The nDim^2 x nNodes matrix, entry @f$ (l\, n_{\text{dim}} + m, B) @f$.
     */
    template < int nDim >
    Eigen::MatrixXd secondGradientsOfNodes( const std::array< Eigen::MatrixXd, nDim >& d2N_dYdY, int nNodes )
    {
      Eigen::MatrixXd secondGradients( nDim * nDim, nNodes );
      for ( int B = 0; B < nNodes; B++ ) {
        const Fastor::Tensor< double, nDim, nDim > d2N_B_dYdY = secondGradientOfNode< nDim >( d2N_dYdY, B );
        for ( int l = 0; l < nDim; l++ )
          for ( int m = 0; m < nDim; m++ )
            secondGradients( l * nDim + m, B ) = d2N_B_dYdY( l, m );
      }
      return secondGradients;
    }

    /**
     * @brief The moment-weighted second gradient @f$ \boldsymbol{D} @f$ of a test function.
     * @tparam nDim Spatial dimension.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle, @f$ \partial Y / \partial x @f$.
     * @param[in] d2T_dYdY The second gradient of the test function of the node.
     * @param[in] secondMoments The second moments @f$ \boldsymbol{M}^Y @f$ of the particle domain.
     * @param[in] J_Y The volume ratio @f$ J_Y @f$ of the intermediate reference configuration.
     * @return @f$ D_{iK} @f$.
     */
    template < int nDim >
    Fastor::Tensor< double, nDim, nDim > momentWeightedSecondGradient(
      const Fastor::Tensor< double, nDim, nDim >& dY_dx,
      const Fastor::Tensor< double, nDim, nDim >& d2T_dYdY,
      const Fastor::Tensor< double, nDim, nDim >& secondMoments,
      double                                      J_Y )
    {
      using namespace Marmot::Math::Tensors::FastorIndices;
      const Fastor::Tensor< double, nDim, nDim > d2T_dYdY_x_M = Fastor::einsum< ij, jk >( d2T_dYdY, secondMoments ) /
                                                                J_Y;
      return Fastor::einsum< ji, jk >( dY_dx, d2T_dYdY_x_M );
    }

    /**
     * @brief Add the derivatives of the stabilization of node @f$ A @f$ with respect to the second gradient of the
     * displacement increment,
     * @f$ \partial f_{Aj} / \partial (\partial^2 \Delta u_k / \partial Y_l \partial Y_m) = D_{im}\,
     * \partial \tau_{ij} / \partial \Delta F_{kl} @f$.
     * @tparam nDim        Spatial dimension.
     * @tparam nDofPerNode Number of degrees of freedom per node; the displacement comes first.
     * @param[in] D The moment-weighted second gradient of the test function of the node.
     * @param[in] dTau_dDeltaF The derivative of the Kirchhoff stress with respect to the increment of the deformation
     * gradient.
     * @param[in,out] rows The derivatives of the node; the stabilization is added.
     */
    template < int nDim, int nDofPerNode >
    void addSecondGradientRows( const Fastor::Tensor< double, nDim, nDim >&             D,
                                const Fastor::Tensor< double, nDim, nDim, nDim, nDim >& dTau_dDeltaF,
                                TangentRowsOfNode< nDim, nDofPerNode >&                 rows )
    {
      for ( int l = 0; l < nDim; l++ )
        for ( int m = 0; m < nDim; m++ )
          for ( int j = 0; j < nDim; j++ )
            for ( int k = 0; k < nDim; k++ ) {
              double dfStabilization = 0.0;
              for ( int i = 0; i < nDim; i++ )
                dfStabilization += D( i, m ) * dTau_dDeltaF( i, j, k, l );
              rows.df_dIncrementSecondGradient[l * nDim + m]( j, k ) += dfStabilization;
            }
    }

    /**
     * @brief Add the derivatives of the stabilization of node @f$ A @f$ through the spatial gradient of the shape
     * functions in the stress gradient of the increment, @f$ \partial f_{Aj} / \partial (\partial \Delta u_k /
     * \partial Y_l) \mathrel{-}= D_{km}\, \partial \tau_{ij} / \partial Y_m\, \partial Y_l / \partial x_i @f$; zero
     * at a zero increment.
     * @tparam nDim        Spatial dimension.
     * @tparam nDofPerNode Number of degrees of freedom per node; the displacement comes first.
     * @param[in] D The moment-weighted second gradient of the test function of the node.
     * @param[in] dTau_dY The gradient of the Kirchhoff stress of the increment along @f$ \boldsymbol{Y} @f$.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle.
     * @param[in,out] rows The derivatives of the node; the stabilization is added.
     */
    template < int nDim, int nDofPerNode >
    void addStressGradientRows( const Fastor::Tensor< double, nDim, nDim >&       D,
                                const Fastor::Tensor< double, nDim, nDim, nDim >& dTau_dY,
                                const Fastor::Tensor< double, nDim, nDim >&       dY_dx,
                                TangentRowsOfNode< nDim, nDofPerNode >&           rows )
    {
      for ( int l = 0; l < nDim; l++ )
        for ( int j = 0; j < nDim; j++ )
          for ( int k = 0; k < nDim; k++ ) {
            double dfStabilization = 0.0;
            for ( int i = 0; i < nDim; i++ )
              for ( int m = 0; m < nDim; m++ )
                dfStabilization += D( k, m ) * dTau_dY( i, j, m ) * dY_dx( l, i );
            rows.df_dIncrementGradient[l]( j, k ) -= dfStabilization;
          }
    }

  } // namespace NSNI

} // namespace Marmot::Discretizations::Meshfree
