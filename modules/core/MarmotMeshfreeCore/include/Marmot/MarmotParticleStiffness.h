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
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotUtils.h"
#include <Eigen/Core>
#include <Fastor/Fastor.h>
#include <algorithm>
#include <array>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief Sees the stiffness of a particle node by node, as the particle describes it
   * (Particle::describeStiffness()).
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
   * By the chain rule, the block of the stiffness that couples the force @f$ \boldsymbol{f}_A @f$ of node @f$ A @f$
   * with the increment of node @f$ B @f$ is
   * @f[
   *   \boldsymbol{K}_{AB} = \frac{\partial \boldsymbol{f}_A}{\partial \Delta\boldsymbol{q}}\, N_B
   *   + \sum_l \frac{\partial \boldsymbol{f}_A}{\partial (\partial \Delta\boldsymbol{q} / \partial Y_l)}\,
   *   \frac{\partial N_B}{\partial Y_l}
   *   + \sum_{l,m} \frac{\partial \boldsymbol{f}_A}{\partial (\partial^2 \Delta\boldsymbol{q} / \partial Y_l
   *   \partial Y_m)}\, \frac{\partial^2 N_B}{\partial Y_l \partial Y_m} .
   * @f]
   * The particle describes the stiffness in this factored form: first the shape functions at its evaluation point
   * (visitShapeFunctions()), once, then for every node @f$ A @f$ the three derivatives of its force, which do not
   * depend on @f$ B @f$ (visitStiffnessOfNode()). A visitor does with them what it needs, e.g. bound the absolute
   * value of the stiffness (AbsoluteStiffnessBound) without forming its @f$ n^2 @f$ blocks, or form them
   * (StiffnessAssembler; both in MarmotParticleStiffnessVisitors.h, for the host).
   *
   * The visitor is not templated, as Particle is not: the dimension and the number of dofs per node are given with the
   * shape functions, and all arrays are column-major. A particle without second gradients gives nullptr for both
   * second-gradient arrays (of the shape functions and of the forces). With @f$ n_d @f$ = nDim and @f$ n @f$ =
   * nDofPerNode:
   * - @p N: 1 x nNodes, @p dN_dY: @f$ n_d @f$ x nNodes, @p d2N_dYdY: @f$ n_d^2 @f$ x nNodes, entry
   *   @f$ (l\, n_d + m, B) @f$ is @f$ \partial^2 N_B / \partial Y_l \partial Y_m @f$;
   * - each derivative of the force of a node is a block of @f$ n @f$ x @f$ n @f$: row @f$ i @f$ is the component of
   *   the force of the node, column @f$ k @f$ the component of the interpolated increment (the fields of a node, in
   *   its dof order). The @f$ n_d @f$ gradient blocks follow each other, block @f$ l @f$ at offset @f$ l\, n^2 @f$,
   *   and so do the @f$ n_d^2 @f$ second-gradient blocks, block @f$ l\, n_d + m @f$ at offset
   *   @f$ (l\, n_d + m)\, n^2 @f$.
   *
   * The arrays are valid during the call only.
   */
  class StiffnessVisitor {
  public:
    virtual ~StiffnessVisitor() = default;

    /**
     * @brief The shape functions at the evaluation point of the particle, before its nodes.
     * @param[in] nDim The spatial dimension.
     * @param[in] nDofPerNode The number of dofs per node.
     * @param[in] nNodes The number of nodes.
     * @param[in] N The trial shape functions @f$ N_B @f$, 1 x nNodes.
     * @param[in] dN_dY Their gradients @f$ \partial N_B / \partial Y_l @f$, nDim x nNodes.
     * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes; nullptr for a particle without second gradients.
     */
    virtual void visitShapeFunctions( int           nDim,
                                      int           nDofPerNode,
                                      int           nNodes,
                                      const double* N,
                                      const double* dN_dY,
                                      const double* d2N_dYdY ) = 0;

    /**
     * @brief The derivatives of the force of node @p A with respect to the interpolated increment.
     * @param[in] A The node.
     * @param[in] df_dIncrement @f$ \partial \boldsymbol{f}_A / \partial \Delta\boldsymbol{q} @f$, one block.
     * @param[in] df_dIncrementGradient @f$ \partial \boldsymbol{f}_A / \partial (\partial \Delta\boldsymbol{q} /
     * \partial Y_l) @f$, nDim blocks.
     * @param[in] df_dIncrementSecondGradient @f$ \partial \boldsymbol{f}_A / \partial (\partial^2
     * \Delta\boldsymbol{q} / \partial Y_l \partial Y_m) @f$, nDim^2 blocks; nullptr for a particle without second
     * gradients.
     */
    virtual void visitStiffnessOfNode( int           A,
                                       const double* df_dIncrement,
                                       const double* df_dIncrementGradient,
                                       const double* df_dIncrementSecondGradient ) = 0;
  };

  /**
   * @brief The stiffness of one node of a particle: the derivatives of the force of the node with respect to the
   * increment the particle interpolates (see StiffnessVisitor), in the layout StiffnessVisitor documents.
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   */
  template < int nDim, int nDofPerNode >
  struct StiffnessOfNode {
    /// @brief A derivative of the force of the node with respect to a component of the interpolated increment.
    using Block = Eigen::Matrix< double, nDofPerNode, nDofPerNode >;

    /// @f$ \partial \boldsymbol{f}_A / \partial \Delta\boldsymbol{q} @f$.
    Block df_dIncrement;
    /// The blocks @f$ \partial \boldsymbol{f}_A / \partial (\partial \Delta\boldsymbol{q} / \partial Y_l) @f$
    /// side by side, block @f$ l @f$ in the columns @f$ l\, n_{\text{dof/node}} @f$ on (see gradientBlock()).
    Eigen::Matrix< double, nDofPerNode, nDim * nDofPerNode > df_dIncrementGradient;
    /// The blocks @f$ \partial \boldsymbol{f}_A / \partial (\partial^2 \Delta\boldsymbol{q} / \partial Y_l
    /// \partial Y_m) @f$ side by side, block @f$ l\, n_{\text{dim}} + m @f$ (see secondGradientBlock()); zero for a
    /// particle without second gradients.
    Eigen::Matrix< double, nDofPerNode, nDim * nDim * nDofPerNode > df_dIncrementSecondGradient;

    /// @brief Set all derivatives to zero.
    void setZero()
    {
      df_dIncrement.setZero();
      df_dIncrementGradient.setZero();
      df_dIncrementSecondGradient.setZero();
    }

    /**
     * @brief The derivative with respect to the gradient of the interpolated increment along @f$ Y_l @f$.
     * @param[in] l The direction.
     * @return Block @f$ l @f$ of df_dIncrementGradient.
     */
    auto gradientBlock( int l ) { return df_dIncrementGradient.template middleCols< nDofPerNode >( nDofPerNode * l ); }

    /**
     * @brief The derivative with respect to the gradient of the interpolated increment along @f$ Y_l @f$.
     * @param[in] l The direction.
     * @return Block @f$ l @f$ of df_dIncrementGradient.
     */
    auto gradientBlock( int l ) const
    {
      return df_dIncrementGradient.template middleCols< nDofPerNode >( nDofPerNode * l );
    }

    /**
     * @brief The derivative with respect to the second gradient of the interpolated increment along
     * @f$ Y_l, Y_m @f$.
     * @param[in] lm The index @f$ l\, n_{\text{dim}} + m @f$.
     * @return Block @p lm of df_dIncrementSecondGradient.
     */
    auto secondGradientBlock( int lm )
    {
      return df_dIncrementSecondGradient.template middleCols< nDofPerNode >( nDofPerNode * lm );
    }

    /**
     * @brief The derivative with respect to the second gradient of the interpolated increment along
     * @f$ Y_l, Y_m @f$.
     * @param[in] lm The index @f$ l\, n_{\text{dim}} + m @f$.
     * @return Block @p lm of df_dIncrementSecondGradient.
     */
    auto secondGradientBlock( int lm ) const
    {
      return df_dIncrementSecondGradient.template middleCols< nDofPerNode >( nDofPerNode * lm );
    }
  };

  /**
   * @brief Add the rows of node @p A to a matrix of a particle: @f$ \boldsymbol{K}_{AB} @f$ for all nodes @f$ B @f$,
   * from the derivatives of the force of the node (StiffnessOfNode) and the trial shape functions at the evaluation
   * point. computePhysicsKernels() of a particle forms its tangent this way, from the stiffness of each node and its
   * inertia.
   *
   * @tparam nDim        Spatial dimension.
   * @tparam nDofPerNode Number of degrees of freedom per node.
   * @param[in] ofNode The derivatives of the force of node @p A with respect to the interpolated increment (in
   * computePhysicsKernels() with its inertia).
   * @param[in] A The node.
   * @param[in] N The trial shape functions @f$ N_B @f$, 1 x nNodes.
   * @param[in] dN_dY Their gradients @f$ \partial N_B / \partial Y_l @f$, nDim x nNodes.
   * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes, entry @f$ (l\, n_{\text{dim}} + m, B) @f$ is
   * @f$ \partial^2 N_B / \partial Y_l \partial Y_m @f$; empty for a particle without second gradients.
   * @param[in,out] K The matrix, column-major, (nNodes nDofPerNode) x (nNodes nDofPerNode); the rows are added.
   */
  template < int nDim, int nDofPerNode >
  void addRowsOfNode( const StiffnessOfNode< nDim, nDofPerNode >& ofNode,
                      int                                         A,
                      const Eigen::MatrixXd&                      N,
                      const Eigen::MatrixXd&                      dN_dY,
                      const Eigen::MatrixXd&                      d2N_dYdY,
                      Eigen::Ref< Eigen::MatrixXd >               K )
  {
    const int idxA = nDofPerNode * A;
    for ( int B = 0; B < N.size(); B++ ) {
      typename StiffnessOfNode< nDim, nDofPerNode >::Block K_AB = ofNode.df_dIncrement * N( B );
      for ( int l = 0; l < nDim; l++ )
        K_AB += ofNode.gradientBlock( l ) * dN_dY( l, B );
      if ( d2N_dYdY.size() )
        for ( int lm = 0; lm < nDim * nDim; lm++ )
          K_AB += ofNode.secondGradientBlock( lm ) * d2N_dYdY( lm, B );
      K.template block< nDofPerNode, nDofPerNode >( idxA, nDofPerNode * B ) += K_AB;
    }
  }

  /**
   * @brief Describe the stiffness of a particle node by node to a visitor: the shape functions, then the stiffness
   * of every node. What the particles' describeStiffness() share, once their material point is evaluated at the state
   * the stiffness is for.
   *
   * @tparam nDim              Spatial dimension.
   * @tparam nDofPerNode       Number of degrees of freedom per node.
   * @tparam AddStiffnessOfNode A callable `(int A, StiffnessOfNode& ofNode)` that adds the stiffness of node @p A to
   * the (zeroed) @p ofNode, e.g. a particle's addStiffnessOfNode().
   * @param[in,out] visitor The visitor.
   * @param[in] N The trial shape functions, 1 x nNodes.
   * @param[in] dN_dY Their gradients, nDim x nNodes.
   * @param[in] d2N_dYdY Their second gradients, nDim^2 x nNodes; empty for a particle without second gradients.
   * @param[in] addStiffnessOfNode The stiffness of a node.
   */
  template < int nDim, int nDofPerNode, typename AddStiffnessOfNode >
  void describeStiffnessNodeByNode( StiffnessVisitor&         visitor,
                                    const Eigen::MatrixXd&    N,
                                    const Eigen::MatrixXd&    dN_dY,
                                    const Eigen::MatrixXd&    d2N_dYdY,
                                    const AddStiffnessOfNode& addStiffnessOfNode )
  {
    const int nNodes = static_cast< int >( N.size() );
    visitor.visitShapeFunctions( nDim,
                                 nDofPerNode,
                                 nNodes,
                                 N.data(),
                                 dN_dY.data(),
                                 d2N_dYdY.size() ? d2N_dYdY.data() : nullptr );
    const bool                           withSecondGradients = d2N_dYdY.size() != 0;
    StiffnessOfNode< nDim, nDofPerNode > ofNode;
    for ( int A = 0; A < nNodes; A++ ) {
      ofNode.setZero();
      addStiffnessOfNode( A, ofNode );
      visitor.visitStiffnessOfNode( A,
                                    ofNode.df_dIncrement.data(),
                                    ofNode.df_dIncrementGradient.data(),
                                    withSecondGradients ? ofNode.df_dIncrementSecondGradient.data() : nullptr );
    }
  }

  /**
   * @brief What an evaluation of a material point changes, saved, and put back when this goes out of scope: the
   * response and the tangents the material point keeps, and its state variables.
   *
   * @details A particle describes its stiffness (Particle::describeStiffness()) at a zero increment of its material
   * point, which overwrites all three. It saves them first; they are put back also if the evaluation or the visitor
   * throws, so the particle is left as it was found. The state variables are those the material point was assigned
   * (MaterialPoint::getAssignedStateVars()). The density the material point keeps is not touched: it is taken once,
   * at initializeYourself().
   *
   * @tparam MaterialPoint A material point with the public members @c response and @c tangents.
   */
  template < typename MaterialPoint >
  class SavedMaterialPoint {
  public:
    /**
     * @brief Save the response, the tangents and the state variables of a material point.
     * @param[in,out] materialPoint The material point; it is put back on destruction.
     */
    explicit SavedMaterialPoint( MaterialPoint& materialPoint )
      : _materialPoint( materialPoint ),
        _response( materialPoint.response ),
        _tangents( materialPoint.tangents ),
        _stateVars( materialPoint.getAssignedStateVars() ),
        _savedStateVars( _stateVars.stateLocation, _stateVars.stateLocation + _stateVars.stateSize )
    {
    }

    /// @brief Put the response, the tangents and the state variables back.
    ~SavedMaterialPoint()
    {
      _materialPoint.response = _response;
      _materialPoint.tangents = _tangents;
      std::copy( _savedStateVars.begin(), _savedStateVars.end(), _stateVars.stateLocation );
    }

    SavedMaterialPoint( const SavedMaterialPoint& )            = delete; ///< not copyable: it puts back once
    SavedMaterialPoint& operator=( const SavedMaterialPoint& ) = delete; ///< not copyable: it puts back once

  private:
    MaterialPoint&                            _materialPoint;  ///< the material point
    const decltype( MaterialPoint::response ) _response;       ///< its saved response
    const decltype( MaterialPoint::tangents ) _tangents;       ///< its saved tangents
    const Marmot::Utilities::StateView        _stateVars;      ///< where its state variables are
    const std::vector< double >               _savedStateVars; ///< their saved values
  };

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
     * @brief The symmetrized smoothed second gradients of all nodes, in the layout of addRowsOfNode().
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
     * @param[in,out] ofNode The derivatives of the node; the stabilization is added.
     */
    template < int nDim, int nDofPerNode >
    void addSecondGradientRows( const Fastor::Tensor< double, nDim, nDim >&             D,
                                const Fastor::Tensor< double, nDim, nDim, nDim, nDim >& dTau_dDeltaF,
                                StiffnessOfNode< nDim, nDofPerNode >&                   ofNode )
    {
      for ( int l = 0; l < nDim; l++ )
        for ( int m = 0; m < nDim; m++ )
          for ( int j = 0; j < nDim; j++ )
            for ( int k = 0; k < nDim; k++ ) {
              double dfStabilization = 0.0;
              for ( int i = 0; i < nDim; i++ )
                dfStabilization += D( i, m ) * dTau_dDeltaF( i, j, k, l );
              ofNode.secondGradientBlock( l * nDim + m )( j, k ) += dfStabilization;
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
     * @param[in,out] ofNode The derivatives of the node; the stabilization is added.
     */
    template < int nDim, int nDofPerNode >
    void addStressGradientRows( const Fastor::Tensor< double, nDim, nDim >&       D,
                                const Fastor::Tensor< double, nDim, nDim, nDim >& dTau_dY,
                                const Fastor::Tensor< double, nDim, nDim >&       dY_dx,
                                StiffnessOfNode< nDim, nDofPerNode >&             ofNode )
    {
      for ( int l = 0; l < nDim; l++ )
        for ( int j = 0; j < nDim; j++ )
          for ( int k = 0; k < nDim; k++ ) {
            double dfStabilization = 0.0;
            for ( int i = 0; i < nDim; i++ )
              for ( int m = 0; m < nDim; m++ )
                dfStabilization += D( k, m ) * dTau_dY( i, j, m ) * dY_dx( l, i );
            ofNode.gradientBlock( l )( j, k ) -= dfStabilization;
          }
    }

  } // namespace NSNI

} // namespace Marmot::Discretizations::Meshfree
