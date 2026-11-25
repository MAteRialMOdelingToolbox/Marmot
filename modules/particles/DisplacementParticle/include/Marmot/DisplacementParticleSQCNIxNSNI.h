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

#include "Marmot/DisplacementParticleSQCNI.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/NewmarkBetaIntegrator.h"
/* Fastor/Fastor.h ends with `#pragma warning( default : ... )`; scope it so it neither overrides the build's warning
 * flags nor leaks into consumers. */
#ifdef _MSC_VER
#  pragma warning( push )
#endif
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#  pragma warning( pop )
#endif

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::DisplacementParticleSQCNIxNSNI
   * @brief Implements a Displacement Particle with Stabilized Quasi-Conforming Nodal Integration (SQCNI)
   *        and Naturally Stabilized Nodal Integration (NSNI) for enhanced stability and accuracy.
   *
   * @details The shape functions and their smoothed gradients are those of DisplacementParticleSQCNI (the smoothing
   * domain follows the SmoothingDomainUpdateType, i.e. this class also provides the SNNI x NSNI variants). In
   * addition, smoothed second derivatives are computed from the gradients at the face centers of the smoothing
   * domain,
   * @f[
   *   N_{B,IJ} = \frac{1}{2}\left( \hat N_{B,IJ} + \hat N_{B,JI} \right), \qquad
   *   \hat N_{B,IJ} = \frac{1}{V_s} \sum_f \frac{\partial N_B}{\partial Y_I}(\boldsymbol{Y}_f)\,(N_J\,dA)_f .
   * @f]
   * The residual of DisplacementParticle is augmented by the NSNI stabilization, the second-order term of a Taylor
   * expansion of the test function gradient and of the stress about the particle center,
   * @f[
   *   r^{\mathrm{stab}}_{Aj} = \Delta F^{-1}_{Ii}\,N_{A,IJ}\,\frac{M_{JK}}{J_Y}\,
   *   \frac{\partial \tau_{ij}}{\partial Y_K}, \qquad
   *   \frac{\partial \tau_{ij}}{\partial Y_K} \approx \frac{\partial \tau_{ij}}{\partial \Delta F_{mM}}\,
   *   \Delta q_{Bm}\,N_{B,MK},
   * @f]
   * with the second moments @f$ M_{JK} = \int_{\Omega_Y} (Y_J - Y_{c,J})(Y_K - Y_{c,K})\,dV @f$ of the particle
   * geometry in the intermediate configuration and @f$ J_Y = \det\boldsymbol{F}_n @f$ (so that
   * @f$ M/J_Y @f$ is referred to the undeformed volume, like @f$ V_0 @f$ in the base residual).
   *
   * The tangent of the stabilization differentiates the second gradient of the displacement increment and the
   * inverse of @f$ \Delta\boldsymbol{F} @f$, but not @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$
   * itself: the term with @f$ \partial^2\boldsymbol{\tau}/\partial\Delta\boldsymbol{F}^2 @f$ is omitted, since the
   * material interface does not provide it. The tangent is therefore approximate by construction (the module test
   * measures a relative error of about @f$ 10^{-3} @f$); all other parts are exact.
   *
   * The second moments are initialized at construction (the intermediate configuration coincides with the undeformed
   * one until the first accepted increment) and updated in acceptStateAndPosition().
   *
   * @tparam nDim The number of spatial dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   */
  template < int nDim, int nVertices >
  class DisplacementParticleSQCNIxNSNI : public DisplacementParticleSQCNI< nDim, nVertices > {

    using ParentPointParticle = DisplacementParticle< nDim >;                 ///< the point particle base class
    using ParentSQCNIParticle = DisplacementParticleSQCNI< nDim, nVertices >; ///< the SQCNI base class

    using TensorD    = Fastor::Tensor< double, nDim >;                        ///< vector of size nDim
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;                  ///< second-order tensor of size nDim
    using TensorDDD  = Fastor::Tensor< double, nDim, nDim, nDim >;            ///< third-order tensor of size nDim
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >;      ///< fourth-order tensor of size nDim

    /**
     * @brief Moments of inertia of the particle in the intermediate reference configuration.
     *
     * This tensor stores the second moments of area/volume of the deformed particle domain
     * which are crucial for the stabilization terms in the NSNI formulation.
     */
    TensorDD _momentsOfInertia_IntermediateReference;

    /// second gradient of the displacement increment, @f$ \partial \Delta F_{iJ}/\partial Y_K @f$ (explicit scheme
    /// only)
    TensorDDD _d2x_dYdY;

    /// velocity gradient with respect to @f$ \boldsymbol{Y} @f$ (explicit scheme only, used by computeLumpedMomentum())
    TensorDD _dv_dY = TensorDD( 0.0 );

    /**
     * @brief Second derivatives of the shape functions with respect to the intermediate coordinates.
     *
     * This array stores the second derivatives of the shape functions, `d^2N / dY_i dY_j`,
     * where `Y` refers to the intermediate coordinates. Each element `_d2N_dYdY[i]`
     * represents `d^2N / dY_i dY`, which is a matrix of size `nDim x nNodes`.
     */
    std::array< Eigen::MatrixXd, nDim > _d2N_dYdY;

    /**
     * @brief The moment-weighted second gradient @f$ \boldsymbol{D} @f$ of the test function of node @p A
     * (NSNI::momentWeightedSecondGradient()): the stabilization force of the node is @f$ D_{iK}\, \partial \tau_{ij}
     * / \partial Y_K @f$.
     * @param[in] A The node.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle, @f$ \partial Y / \partial x @f$.
     * @return @f$ \boldsymbol{D} @f$.
     */
    TensorDD momentWeightedSecondGradientOfTestFunction( int A, const TensorDD& dY_dx ) const
    {
      return NSNI::momentWeightedSecondGradient< nDim >( dY_dx,
                                                         NSNI::secondGradientOfNode< nDim >( _d2N_dYdY, A ),
                                                         _momentsOfInertia_IntermediateReference,
                                                         Fastor::determinant( this->_mp->dY_dX() ) );
    }

    /**
     * @brief The stiffness of node @p A: those of DisplacementParticle plus the stabilization at the current
     * state, which depends on the second gradient of the increment,
     * @f[
     *   \frac{\partial f_{Aj}}{\partial (\partial^2 \Delta u_k / \partial Y_l \partial Y_m)} = D_{im}
     *   \frac{\partial \tau_{ij}}{\partial \Delta F_{kl}}
     * @f]
     * (see momentWeightedSecondGradientOfTestFunction()). The stabilization's dependence on the spatial gradient of
     * the shape functions through the stress gradient of the increment vanishes at a zero increment; it is added in
     * computePhysicsKernels().
     * @param[in] A The node.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle, @f$ \partial Y / \partial x @f$,
     * the same for all its nodes: computed once per particle by the caller.
     * @param[in,out] ofNode The derivatives; the stiffness of the node is added.
     */
    void addStiffnessOfNode( int                                         A,
                             const Fastor::Tensor< double, nDim, nDim >& dY_dx,
                             StiffnessOfNode< nDim, nDim >&              ofNode ) const override;

    /**
     * @brief The smoothed second gradients of the trial shape functions, symmetrized, in the layout of
     * addRowsOfNode() (NSNI::secondGradientsOfNodes()).
     * @return The nDim^2 x nNodes matrix.
     */
    Eigen::MatrixXd secondGradientsOfTrialFunctions() const override
    {
      return NSNI::secondGradientsOfNodes< nDim >( _d2N_dYdY, this->_nNodes );
    }

  public:
    /**
     * @brief Constructor for DisplacementParticleSQCNIxNSNI.
     *
     * @param[in] elementID Unique identifier for the particle.
     * @param[in] nodeCoordinates Pointer to an array of node coordinates defining the particle's initial geometry.
     * @param[in] nNodeCoordiantes Number of node coordinates.
     * @param[in] volume Initial volume of the particle; must be 0 (see DisplacementParticleSQCNI).
     * @param[in] materialName Name of the material assigned to the particle.
     * @param[in] materialProperties Pointer to an array of material properties.
     * @param[in] sizeMaterialProperties Size of the material properties array.
     * @param[in] approximation Reference to the meshfree approximation object.
     * @param[in] smoothingVolumeUpdateType Type of update strategy for the smoothing domain volume.
     * @throws std::invalid_argument if the provided volume is not zero.
     */
    DisplacementParticleSQCNIxNSNI( int                                                     elementID,
                                    const double*                                           nodeCoordinates,
                                    int                                                     nNodeCoordiantes,
                                    double                                                  volume,
                                    const std::string&                                      materialName,
                                    const double*                                           materialProperties,
                                    int                                                     sizeMaterialProperties,
                                    const Marmot::Discretizations::Meshfree::Approximation& approximation,
                                    ParentSQCNIParticle::SmoothingDomainUpdateType          smoothingVolumeUpdateType );

    /**
     * @brief Evaluate the shape functions with the assigned kernel functions, with smoothed first and second
     * derivatives.
     *
     * This method overrides the base class implementation to compute not only the
     * shape functions and their first derivatives but also the second derivatives
     * (`_d2N_dYdY`) required for the NSNI formulation. It uses a boundary integral
     * approach to compute the derivatives.
     */
    void updateShapeFunctions() override
    {
      Eigen::Matrix< double, nDim, 1 > centerCoordinates;
      this->getCenterCoordinates( centerCoordinates.data() );

      ParentPointParticle::_N = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( centerCoordinates.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         ParentPointParticle::_N.data() );

      ParentSQCNIParticle::SectionIntegrals::
        smoothShapeFunctionGradients( this->smoothingDomain(),
                                      ParentPointParticle::_meshfreeApproximation,
                                      ParentPointParticle::_assignedKernelFunctions,
                                      ParentPointParticle::_N,
                                      ParentPointParticle::_dN_dY,
                                      &_d2N_dYdY );

      ParentPointParticle::_T     = ParentPointParticle::_N;
      ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
    }

    /**
     * @brief Computes the internal force vector and tangent stiffness matrix for the particle.
     *
     * This method implements the core physics computations for the NSNI particle.
     * It calculates the internal forces (`fInt`) and the tangent stiffness matrix (`dFInt_ddQ`)
     * based on the current deformation, material response, and the NSNI stabilization terms (see the class
     * description; the tangent of the stabilization is approximate). Inertia is included as in DisplacementParticle.
     *
     * @param[in] dQ Pointer to the incremental nodal displacement vector.
     * @param[in,out] fInt Pointer to the output internal force vector (the contribution is added).
     * @param[in,out] dFInt_ddQ Pointer to the output tangent stiffness matrix (the contribution is added).
     * @param[in] timeNew Current simulation time.
     * @param[in] dT Time step size.
     */
    void computePhysicsKernels( const double* dQ, double* fInt, double* dFInt_ddQ, double timeNew, double dT ) override;

    /**
     * @brief Explicit scheme: updates the material point with the increment dQ, stores the second gradient of the
     * increment, and sets @f$ \boldsymbol{v}_{n+1} = \Delta\boldsymbol{u}/\Delta t @f$,
     * @f$ \boldsymbol{a} = (\boldsymbol{v}_{n+1} - \boldsymbol{v}_n)/\Delta t @f$ and the velocity gradient
     * (skipped for @f$ \Delta t \le 10^{-16} @f$).
     * @param[in] dQ Pointer to the incremental nodal displacement vector.
     * @param[in] timeNew Current simulation time.
     * @param[in] dT Time step size.
     */
    void updatePhysicsExplicit( const double* dQ, double timeNew, double dT ) override;

    /**
     * @brief Explicit scheme: internal force of the current state (as computePhysicsKernels(), without inertia and
     * tangent), including the NSNI stabilization.
     * @param[in,out] fInt Pointer to the internal force vector (the contribution is added).
     */
    void computePhysicsKernelsExplicit( double* fInt ) override;

    /**
     * @brief Lumped momentum, @f$ p_{Ai} \mathrel{+}= \rho_0\,T_A\,V_0\,v_i
     * + \rho_0\,(\partial v_i/\partial Y_J)\,M_{JK}\,(\partial T_A/\partial Y_K)/J_Y @f$, with the velocity gradient
     * of the last updatePhysicsExplicit().
     * @param[in,out] pLumped Lumped momentum vector (nDim values per node), the contribution is added.
     */
    void computeLumpedMomentum( double* pLumped ) const override;

    /**
     * @brief Accepts the current state and position, updating internal variables.
     *
     * This method is called to finalize the state after a successful time step.
     * It updates the base class state and also computes the moments of inertia
     * for the intermediate reference configuration, which are used in the NSNI stabilization.
     */
    void acceptStateAndPosition() override
    {
      ParentSQCNIParticle::acceptStateAndPosition();
      _momentsOfInertia_IntermediateReference = TensorDD( this->_particleDomain.getGeometrySecondMoments().data() );
    };
  };

  template < int nDim, int nVertices >
  DisplacementParticleSQCNIxNSNI< nDim, nVertices >::DisplacementParticleSQCNIxNSNI(
    int                                                     elementID,
    const double*                                           vertexCoordinates,
    int                                                     nVertexCoordinates,
    double                                                  volume,
    const std::string&                                      materialName,
    const double*                                           materialProperties,
    int                                                     sizeMaterialProperties,
    const Marmot::Discretizations::Meshfree::Approximation& approximation,
    ParentSQCNIParticle::SmoothingDomainUpdateType          smoothingVolumeUpdateType )
    : DisplacementParticleSQCNI< nDim, nVertices >( elementID,
                                                    vertexCoordinates,
                                                    nVertexCoordinates,
                                                    volume,
                                                    materialName,
                                                    materialProperties,
                                                    sizeMaterialProperties,
                                                    approximation,
                                                    smoothingVolumeUpdateType )
  {
    // the intermediate reference coincides with the undeformed one until the first accepted increment;
    // without this, the stabilization of the first increment is scaled by uninitialized memory
    _momentsOfInertia_IntermediateReference = TensorDD( this->_particleDomain.getGeometrySecondMoments().data() );
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxNSNI< nDim, nVertices >::computePhysicsKernels( const double* dQ,
                                                                                 double*       fInt,
                                                                                 double*       dFInt_ddQ,
                                                                                 double        timeNew,
                                                                                 double        dT )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;
    using ijmM = Index< i_, j_, m_, M_ >;
    using mMK  = Index< m_, M_, K_ >;
    using ijK  = Fastor::Index< i_, j_, K_ >;

    const auto& _nNodes = this->_nNodes;
    const auto& _N      = this->_N;
    const auto& _dN_dY  = this->_dN_dY;
    const auto& _T      = this->_T;
    const auto& _dT_dY  = this->_dT_dY;
    auto&       _mp     = this->_mp;

    const static TensorDD I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    TensorD   du( 0.0 );
    TensorDD  du_dY( 0.0 );
    TensorDDD d2x_dYdY( 0.0 );

    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;

      const double N_B        = _N( B );
      const auto   dN_B_dY    = TensorD( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen
      const auto   d2N_B_dYdY = NSNI::secondGradientOfNode< nDim >( _d2N_dYdY, B );

      const auto dQU_B = TensorD( dQ + idxB_u );

      du += N_B * dQU_B;

      du_dY += einsum< i, j >( dQU_B, dN_B_dY );

      d2x_dYdY += einsum< i, jk >( dQU_B, d2N_B_dYdY );
    }

    _mp->prepareYourself( timeNew, dT );
    _mp->incrementDeformation( du, du_dY );
    _mp->computeYourself( timeNew, dT );

    const double density0 = _mp->getDensityUndeformed();

    auto v = _mp->getVelocity();
    auto a = _mp->getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::ContinuumMechanics::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                                                 v.data(),
                                                                                 a.data(),
                                                                                 dT,
                                                                                 this->_newmark_beta,
                                                                                 this->_newmark_gamma,
                                                                                 da_ddu.data() );
    _mp->setVelocity( v );
    _mp->setAcceleration( a );

    TensorD r_U( 0.0 );

    const auto& S = _mp->response.S;

    const double V0 = this->getVolumeUndeformed();

    const auto& t = _mp->tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    const auto dY_dx = evaluate( inv( _mp->dx_dY() ) );

    const auto dS_dY = evaluate( einsum< ijmM, mMK >( t.dS_dDeltaF, d2x_dYdY ) );

    const Eigen::MatrixXd d2N_dYdY = secondGradientsOfTrialFunctions();

    StiffnessOfNode< nDim, nDim > rowsOfNode;
    for ( int A = 0; A < _nNodes; A++ ) {

      const double  T_A     = _T( A );
      const auto    dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const TensorD dT_A_dx = einsum< ji, j >( dY_dx, dT_A_dY );

      r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

      // add inertia
      r_U += density0 * a * T_A * V0;

      const TensorDD d2NA_dxdY_x_MOIScaled = momentWeightedSecondGradientOfTestFunction( A, dY_dx );

      TensorD rU_Stab = einsum< iK, ijK >( d2NA_dxdY_x_MOIScaled, dS_dY );

      r_U += rU_Stab;

      P.template segment< nDim >( nodeBlockSize * A ) += Eigen::Map< Eigen::Matrix< double, nDim, 1 > >( r_U.data() );

      // the rows of the tangent: the stiffness of the node plus its inertia and the stabilization of the increment
      rowsOfNode.setZero();
      addStiffnessOfNode( A, dY_dx, rowsOfNode );
      rowsOfNode.df_dIncrement += density0 * T_A * V0 *
                                  Eigen::Map< const Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >(
                                    da_ddu.data() );

      // the stabilization through the stress gradient of the increment: D_kK dTau_ij/dY_K dN_B/dx_i
      NSNI::addStressGradientRows< nDim >( d2NA_dxdY_x_MOIScaled, dS_dY, dY_dx, rowsOfNode );

      addRowsOfNode( rowsOfNode, A, _N, _dN_dY, d2N_dYdY, K );
    }
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxNSNI< nDim, nVertices >::addStiffnessOfNode(
    int                                         A,
    const Fastor::Tensor< double, nDim, nDim >& dY_dx,
    StiffnessOfNode< nDim, nDim >&              ofNode ) const
  {
    ParentPointParticle::addStiffnessOfNode( A, dY_dx, ofNode );
    NSNI::addSecondGradientRows< nDim >( momentWeightedSecondGradientOfTestFunction( A, dY_dx ),
                                         this->_mp->tangents.dS_dDeltaF,
                                         ofNode );
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxNSNI< nDim, nVertices >::computePhysicsKernelsExplicit( double* fInt )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;
    using ijmM = Index< i_, j_, m_, M_ >;
    using mMK  = Index< m_, M_, K_ >;
    using ijK  = Fastor::Index< i_, j_, K_ >;

    const auto& _nNodes = this->_nNodes;
    const auto& _dT_dY  = this->_dT_dY;
    auto&       _mp     = this->_mp;

    const static TensorDD I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    TensorD r_U( 0.0 );

    const auto& S = _mp->response.S;

    const double V0 = this->getVolumeUndeformed();

    const auto& t = _mp->tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );

    const auto dY_dx = evaluate( inv( _mp->dx_dY() ) );

    const auto dS_dY = evaluate( einsum< ijmM, mMK >( t.dS_dDeltaF, _d2x_dYdY ) );

    // clang-format off
    for ( int A = 0; A < _nNodes; A++ ) {

      const auto  dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const TensorD dT_A_dx = einsum< ji, j >( dY_dx , dT_A_dY );


      const int idxA_u = nodeBlockSize * A;

      r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

      const TensorDD d2NA_dxdY_x_MOIScaled = momentWeightedSecondGradientOfTestFunction( A, dY_dx );

      TensorD rU_Stab = einsum< iK, ijK >( d2NA_dxdY_x_MOIScaled, dS_dY );

      r_U += rU_Stab;

      {
        using namespace Eigen;
        P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
      }
    }
    // clang-format on
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxNSNI< nDim, nVertices >::updatePhysicsExplicit( const double* dQ,
                                                                                 double        timeNew,
                                                                                 double        dT )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    const auto& _nNodes = this->_nNodes;
    const auto& _N      = this->_N;
    const auto& _dN_dY  = this->_dN_dY;
    auto&       _mp     = this->_mp;

    const static TensorDD I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    TensorD  du( 0.0 );
    TensorDD du_dY( 0.0 );

    _d2x_dYdY.zeros();
    for ( int B = 0; B < _nNodes; B++ ) {

      const int idxB_u = nodeBlockSize * B;

      const double N_B        = _N( B );
      const auto   dN_B_dY    = TensorD( _dN_dY.col( B ).data() ); // works because ColumnMajor of Eigen
      const auto   d2N_B_dYdY = NSNI::secondGradientOfNode< nDim >( _d2N_dYdY, B );

      const auto dQU_B = TensorD( dQ + idxB_u );

      du += N_B * dQU_B;

      du_dY += einsum< i, j >( dQU_B, dN_B_dY );

      // the second gradient of x = Y + du is that of du: Y is the coordinate itself
      _d2x_dYdY += einsum< i, jk >( dQU_B, d2N_B_dYdY );
    }

    _mp->prepareYourself( timeNew, dT );
    _mp->incrementDeformation( du, du_dY );
    _mp->computeYourself( timeNew, dT );
    // a zero-length increment (the start of a step) leaves the rates as they are
    if ( dT <= 0.0 )
      return;
    const auto    v_n  = _mp->getVelocity();
    const TensorD v_np = du / dT;
    _mp->setVelocity( v_np );
    _mp->setAcceleration( evaluate( v_np - v_n ) / dT );

    _dv_dY = du_dY / dT;
  }

  template < int nDim, int nVertices >
  void DisplacementParticleSQCNIxNSNI< nDim, nVertices >::computeLumpedMomentum( double* pLumped ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    const double density0 = this->_mp->getDensityUndeformed();
    const double V0       = this->getVolumeUndeformed();
    const auto   v        = this->_mp->getVelocity();
    const auto&  _dT_dY   = this->_dT_dY;

    const TensorDD _dv_dY_x_Y2 = einsum< ij, jk >( _dv_dY, _momentsOfInertia_IntermediateReference );
    //
    TensorD x_p;
    this->getCenterCoordinates( x_p.data() );

    for ( int A = 0; A < this->_nNodes; A++ ) { // Use base class _nNodes

      const double T_A     = this->_T( A );     // Use base class _T
      const auto   dT_A_dY = Tensor< double, nDim >( _dT_dY.col( A ).data() );

      const int idxA_u = this->nDofPerNodeU * A;

      const TensorD aux = einsum< ij, j >( _dv_dY_x_Y2, dT_A_dY );

      const double detJIntermediate = determinant( this->_mp->dY_dX() );

      for ( int i = 0; i < this->nDofPerNodeU; i++ ) {
        pLumped[idxA_u + i] += density0 * T_A * V0 * v[i];
        pLumped[idxA_u + i] += density0 * aux[i] / detJIntermediate;
      }
    }
  }

} // namespace Marmot::Discretizations::Meshfree
