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

#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleStiffness.h"
#include <Eigen/Core>
#include <vector>

/**
 * @file
 * @brief The visitors of a particle's stiffness that a host applies (Particle::describeStiffness()): the bound of its
 * absolute value, for the stable time increment of an explicit scheme, and its assembly. The particles do not
 * include this header; both visitors and their free functions are compiled in Marmot, for nDim 2 and 3 (a host that
 * compiled them with its own flags, e.g. with FMA contraction, would change the last bits of their results).
 */

namespace Marmot::Discretizations::Meshfree {

  /// The largest number of dofs per node the visitors of this header are compiled for (a node's rows are read as
  /// blocks of a size known at compile time); a particle with more is refused.
  inline constexpr int maxDofsPerNode = 12;

  /**
   * @brief Bounds the absolute value of a particle's stiffness times nonnegative weights,
   * @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$, from the stiffness it describes node by node, without forming it.
   *
   * @details What a stable time increment needs from a stiffness is a bound of its spectral radius, e.g. the
   * Collatz-Wielandt bound, @f$ \max_i (|\boldsymbol{K}|\boldsymbol{w})_i / w_i @f$ for positive weights
   * @f$ \boldsymbol{w} @f$. The stiffness of a particle with @f$ n @f$ nodes has @f$ n^2 @f$ blocks; this bound of
   * @f$ |\boldsymbol{K}|\boldsymbol{W} @f$ costs @f$ O(n) @f$: the rows of node @f$ A @f$ are linear in the shape
   * functions of the other node @f$ B @f$ (see StiffnessVisitor), and the sum over @f$ B @f$ is taken from a few
   * weighted moments of the shape functions, computed once per particle (visitShapeFunctions()). Row @f$ i @f$ of
   * node @f$ A @f$, column component @f$ k @f$, with the weights @f$ w_B @f$ of component @f$ k @f$ of the nodes and
   * the factors @f$ a = \partial f_{Ai}/\partial \Delta q_k @f$, @f$ \boldsymbol{b} @f$ (its gradient factors) and
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
   * spread evenly. The terms are added by the triangle inequality. The bound holds for every stiffness of the form of
   * StiffnessVisitor, symmetric or not; the second gradients of the shape functions must be symmetric (their
   * symmetric components are used).
   *
   * A host uses it through addBoundOfAbsoluteStiffnessTimes().
   *
   * @tparam nDim Spatial dimension.
   */
  template < int nDim >
  class AbsoluteStiffnessBound : public StiffnessVisitor {
  public:
    /// The number of independent components of a symmetric second gradient.
    static constexpr int nSym = nDim * ( nDim + 1 ) / 2;

    /**
     * @brief A bound to be added to @p bound, for the weights @p columnWeights.
     * @param[in] columnWeights The nonnegative weights @f$ \boldsymbol{W} @f$, column-major, (nNodes nDofPerNode) x
     * nColumns, in the dof order of the particle.
     * @param[in] nColumns The number of columns of @p columnWeights.
     * @param[in,out] bound The bound, column-major, (nNodes nDofPerNode) x nColumns; the bound of every row is added.
     */
    AbsoluteStiffnessBound( const double* columnWeights, int nColumns, double* bound )
      : _columnWeights( columnWeights ), _nColumns( nColumns ), _bound( bound )
    {
    }

    /**
     * @brief Compute the weighted moments of the shape functions.
     * @param[in] nDimOfParticle The spatial dimension of the particle; must be nDim.
     * @param[in] nDofPerNode The number of dofs per node.
     * @param[in] nNodes The number of nodes.
     * @param[in] N The trial shape functions, 1 x nNodes.
     * @param[in] dN_dY Their gradients, nDim x nNodes.
     * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes; nullptr for a particle without second gradients.
     * @throws std::invalid_argument if the particle is not of dimension nDim, or has more than maxDofsPerNode dofs per
     * node.
     */
    void visitShapeFunctions( int           nDimOfParticle,
                              int           nDofPerNode,
                              int           nNodes,
                              const double* N,
                              const double* dN_dY,
                              const double* d2N_dYdY ) override;

    /**
     * @brief Add the bound of the rows of node @p A: @f$ (|\boldsymbol{K}|\boldsymbol{W})_{Ai,c} @f$ for all
     * components @f$ i @f$ and columns @f$ c @f$.
     * @param[in] A The node.
     * @param[in] df_dIncrement The derivative of its force with respect to the interpolated increment.
     * @param[in] df_dIncrementGradient The derivatives with respect to its gradient.
     * @param[in] df_dIncrementSecondGradient The derivatives with respect to its second gradient.
     */
    void visitStiffnessOfNode( int           A,
                               const double* df_dIncrement,
                               const double* df_dIncrementGradient,
                               const double* df_dIncrementSecondGradient ) override;

  private:
    /**
     * @brief Add the bound of the rows of node @p A (see visitStiffnessOfNode()), with the number of dofs per node
     * known at compile time: the derivatives are read as typed blocks, and the loops over the components unrolled.
     * @tparam nDofPerNode The number of dofs per node.
     * @param[in] A The node.
     * @param[in] df_dIncrement The derivative of its force with respect to the interpolated increment.
     * @param[in] df_dIncrementGradient The derivatives with respect to its gradient.
     * @param[in] df_dIncrementSecondGradient The derivatives with respect to its second gradient.
     */
    template < int nDofPerNode >
    void addBoundOfRowsOfNode( int           A,
                               const double* df_dIncrement,
                               const double* df_dIncrementGradient,
                               const double* df_dIncrementSecondGradient );

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
                            Eigen::Matrix< double, n, n >&       moments );

    /**
     * @brief The independent components @f$ (l, m), l \le m @f$ of a symmetric second gradient.
     * @param[in] d2N_dYdY The second gradient of one node, entry @f$ l\, n_{\text{dim}} + m @f$.
     * @return The nSym components, in the order @f$ (0,0), (0,1), \ldots, (1,1), \ldots @f$.
     */
    static Eigen::Matrix< double, nSym, 1 > symmetricComponents( const double* d2N_dYdY );

    const double*          _columnWeights;               ///< the weights, column-major
    const int              _nColumns;                    ///< the number of columns of the weights
    double*                _bound;                       ///< the bound, column-major
    int                    _nDof                = 0;     ///< the number of dofs of the particle
    bool                   _withSecondGradients = false; ///< the particle has second gradients
    std::vector< Moments > _moments;                     ///< entry @f$ c\, n_{\text{dof/node}} + k @f$
    /// @cond INTERNAL (a pointer to a member function, which the Sphinx C++ domain cannot parse)
    /// addBoundOfRowsOfNode() for the dofs per node of the particle, chosen by visitShapeFunctions()
    void ( AbsoluteStiffnessBound::*_addBoundOfRowsOfNode )( int,
                                                             const double*,
                                                             const double*,
                                                             const double* ) = nullptr;
    /// @endcond
  };

  /**
   * @brief Add an upper bound of @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$, the absolute value of a particle's
   * stiffness times nonnegative weights, without forming the stiffness (AbsoluteStiffnessBound, which the particle
   * describes its stiffness to; see Particle::describeStiffness() for the state it is of). The particle is left as it
   * was found.
   *
   * @details The bound a host needs for the stable time increment of an explicit scheme, e.g. by the Collatz-Wielandt
   * bound of the mass-scaled stiffness, which needs @f$ |\boldsymbol{K}|\boldsymbol{w} @f$ only. It costs O(nodes)
   * instead of the O(nodes^2) of the stiffness.
   * @param[in,out] particle The particle (Particle::describeStiffness()).
   * @param[in] columnWeights The nonnegative weights @f$ \boldsymbol{W} @f$, column-major, nDof x nColumns, in the dof
   * order of the particle.
   * @param[in] nColumns The number of columns of @p columnWeights.
   * @param[in,out] bound Column-major, nDof x nColumns: the bound is added.
   * @param[in] timeNew The time of the zero increment the stiffness is evaluated with.
   * @param[in] dT The time increment of the zero increment.
   * @throws std::invalid_argument for a particle of a dimension other than 2 or 3, or with more than maxDofsPerNode
   * dofs per node; std::runtime_error from a particle that does not describe its stiffness.
   */
  void addBoundOfAbsoluteStiffnessTimes( Particle&     particle,
                                         const double* columnWeights,
                                         int           nColumns,
                                         double*       bound,
                                         double        timeNew,
                                         double        dT );

  /**
   * @brief Assembles the stiffness of a particle from the stiffness it describes node by node: the blocks
   * @f$ \boldsymbol{K}_{AB} @f$ of all nodes (see StiffnessVisitor), with addRowsOfNode(), as the particle's
   * computePhysicsKernels() forms its tangent, without the inertia.
   *
   * A host uses it through assembleStiffness().
   *
   * @tparam nDim Spatial dimension.
   */
  template < int nDim >
  class StiffnessAssembler : public StiffnessVisitor {
  public:
    /**
     * @brief A stiffness to be added to @p K.
     * @param[in,out] K The stiffness, column-major, nDof x nDof; the stiffness of the particle is added.
     */
    explicit StiffnessAssembler( double* K ) : _K( K ) {}

    /**
     * @brief Keep the shape functions for the rows of the nodes.
     * @param[in] nDimOfParticle The spatial dimension of the particle; must be nDim.
     * @param[in] nDofPerNode The number of dofs per node.
     * @param[in] nNodes The number of nodes.
     * @param[in] N The trial shape functions, 1 x nNodes.
     * @param[in] dN_dY Their gradients, nDim x nNodes.
     * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes; nullptr for a particle without second gradients.
     * @throws std::invalid_argument if the particle is not of dimension nDim, or has more than maxDofsPerNode dofs per
     * node.
     */
    void visitShapeFunctions( int           nDimOfParticle,
                              int           nDofPerNode,
                              int           nNodes,
                              const double* N,
                              const double* dN_dY,
                              const double* d2N_dYdY ) override;

    /**
     * @brief Add the rows of node @p A: @f$ \boldsymbol{K}_{AB} @f$ for all nodes @f$ B @f$.
     * @param[in] A The node.
     * @param[in] df_dIncrement The derivative of its force with respect to the interpolated increment.
     * @param[in] df_dIncrementGradient The derivatives with respect to its gradient.
     * @param[in] df_dIncrementSecondGradient The derivatives with respect to its second gradient.
     */
    void visitStiffnessOfNode( int           A,
                               const double* df_dIncrement,
                               const double* df_dIncrementGradient,
                               const double* df_dIncrementSecondGradient ) override;

  private:
    /**
     * @brief Add the rows of node @p A (see visitStiffnessOfNode()) with addRowsOfNode(), for the number of dofs per
     * node at compile time.
     * @tparam nDofPerNode The number of dofs per node.
     * @param[in] A The node.
     * @param[in] df_dIncrement The derivative of its force with respect to the interpolated increment.
     * @param[in] df_dIncrementGradient The derivatives with respect to its gradient.
     * @param[in] df_dIncrementSecondGradient The derivatives with respect to its second gradient.
     */
    template < int nDofPerNode >
    void addRowsOfNodeToK( int           A,
                           const double* df_dIncrement,
                           const double* df_dIncrementGradient,
                           const double* df_dIncrementSecondGradient );

    double*         _K;        ///< the stiffness, column-major
    int             _nDof = 0; ///< the number of dofs of the particle
    Eigen::MatrixXd _N;        ///< the trial shape functions, 1 x nNodes
    Eigen::MatrixXd _dN_dY;    ///< their gradients, nDim x nNodes
    Eigen::MatrixXd _d2N_dYdY; ///< their second gradients, nDim^2 x nNodes; empty without second gradients
    /// @cond INTERNAL (a pointer to a member function, which the Sphinx C++ domain cannot parse)
    /// addRowsOfNodeToK() for the dofs per node of the particle, chosen by visitShapeFunctions()
    void ( StiffnessAssembler::*_addRowsOfNodeToK )( int, const double*, const double*, const double* ) = nullptr;
    /// @endcond
  };

  /// @cond INTERNAL (compiled once, in Marmot: src/MarmotParticleStiffnessVisitors.cpp)
  extern template class AbsoluteStiffnessBound< 2 >;
  extern template class AbsoluteStiffnessBound< 3 >;
  extern template class StiffnessAssembler< 2 >;
  extern template class StiffnessAssembler< 3 >;
  /// @endcond

  /**
   * @brief Add the stiffness of a particle to a matrix (StiffnessAssembler, which the particle describes its
   * stiffness to; see Particle::describeStiffness() for the state it is of): the tangent of computePhysicsKernels() at
   * a zero increment, without the inertia. The particle is left as it was found.
   * @param[in,out] particle The particle (Particle::describeStiffness()).
   * @param[in,out] K Column-major, nDof x nDof, in the dof order of the particle: the stiffness is added.
   * @param[in] timeNew The time of the zero increment the stiffness is evaluated with.
   * @param[in] dT The time increment of the zero increment.
   * @throws std::invalid_argument for a particle of a dimension other than 2 or 3, or with more than maxDofsPerNode
   * dofs per node; std::runtime_error from a particle that does not describe its stiffness.
   */
  void assembleStiffness( Particle& particle, double* K, double timeNew, double dT );

} // namespace Marmot::Discretizations::Meshfree
