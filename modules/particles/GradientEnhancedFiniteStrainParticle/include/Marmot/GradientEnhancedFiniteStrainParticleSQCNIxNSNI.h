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

#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNI.h"
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
   * @class Marmot::Discretizations::Meshfree::GradientEnhancedFiniteStrainParticleSQCNIxNSNI
   * @brief SQCNI / SNNI particle with naturally stabilized nodal integration (NSNI) for gradient-enhanced
   *        finite-strain materials.
   *
   * On top of the smoothed gradients of GradientEnhancedFiniteStrainParticleSQCNI, the momentum residual receives the
   * NSNI stabilization term: the integrand is expanded to first order about the particle
   * center, and the second-order term of the integral is kept, with the second moments of the particle domain. The
   * second derivatives of the shape functions are smoothed like the first ones, by the gradients at the face centers,
   * @f[
   *   \frac{\partial^2 N_B}{\partial Y_I\,\partial Y_J} \approx \frac{1}{2V_{\Omega_Y}} \sum_f \Bigl(
   *   \frac{\partial N_B}{\partial Y_I}(\boldsymbol{Y}_f)\,n_J + \frac{\partial N_B}{\partial
   * Y_J}(\boldsymbol{Y}_f)\,n_I \Bigr) dA_f
   * @f]
   * (symmetrized), and the Kirchhoff stress gradient along @f$ \boldsymbol{Y} @f$ is linearized as
   * @f[
   *   \frac{\partial\tau_{ij}}{\partial Y_K} \approx \frac{\partial\tau_{ij}}{\partial\Delta F_{mM}}\,
   *   \frac{\partial^2\Delta x_m}{\partial Y_M\,\partial Y_K}, \qquad
   *   \frac{\partial^2\Delta x_m}{\partial Y_M\,\partial Y_K} = \Delta q^U_{Bm}\,
   *   \frac{\partial^2 N_B}{\partial Y_M\,\partial Y_K} .
   * @f]
   * The stabilization residual added to @f$ r^U_{Aj} @f$ is
   * @f[
   *   r^{U,\mathrm{stab}}_{Aj} = \frac{\partial}{\partial x_i}\Bigl(\frac{\partial N_A}{\partial Y_L}\Bigr)
   *   \frac{M^Y_{LK}}{J_Y}\,\frac{\partial\tau_{ij}}{\partial Y_K}, \qquad
   *   \frac{\partial}{\partial x_i}\Bigl(\frac{\partial N_A}{\partial Y_L}\Bigr) =
   *   \Delta F^{-1}_{mi}\,\frac{\partial^2 N_A}{\partial Y_m\,\partial Y_L},
   * @f]
   * with the second moments @f$ \boldsymbol{M}^Y = J_Y\,\boldsymbol{F}_n\,\boldsymbol{M}^0\,\boldsymbol{F}_n^\mathsf{T}
   * @f$ of the particle domain about its centroid in the intermediate reference configuration (updated at each
   * accepted increment from those, @f$ \boldsymbol{M}^0 @f$, of the undeformed cell) and @f$ J_Y = \det\boldsymbol{F}_n
   * @f$. Note that the stabilization uses the trial second derivatives for both the test and the trial side, i.e.
   * it is not affected by VCI.
   *
   * The stabilization acts on the momentum balance only; the nonlocal residual and its rows of the tangent are those
   * of the SQCNI particle. The stress gradient does not contain the contribution
   * @f$ \partial\boldsymbol{\tau}/\partial\bar{N}\,\nabla_Y\Delta\bar{N} @f$ of the nonlocal field.
   *
   * The tangent of the stabilization term is APPROXIMATE by construction: it differentiates
   * @f$ \partial^2\Delta\boldsymbol{x}/\partial\boldsymbol{Y}^2 @f$ and @f$ \Delta\boldsymbol{F}^{-1} @f$, but not
   * @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$ itself, i.e. it omits
   * @f$ \partial^2\boldsymbol{\tau}/\partial\boldsymbol{F}^2 @f$ (not exposed by the material interface) and the
   * dependence of @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$ on @f$ \bar{N} @f$; the U-N block
   * therefore has no stabilization contribution at all. The nonlocal rows are exact; the module test bounds the
   * error of the U-U block and does not check the U-N block for this particle.
   *
   * @tparam nDim      Spatial dimension (2: plane strain, 3: 3D).
   * @tparam nVertices Number of vertices of the smoothing domain (4: Quad, 8: Hexa).
   * @tparam nNonlocal Number of nonlocal variables (1: a scalar nonlocal damage driver, 6: a nonlocal strain).
   */
  template < int nDim, int nVertices, int nNonlocal >
  class GradientEnhancedFiniteStrainParticleSQCNIxNSNI
    : public GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal > {

    using ParentPointParticle = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >; ///< point particle base
    using ParentSQCNIParticle = GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal >; ///< SQCNI base
    using LagrangeCellType    = ParentSQCNIParticle::LagrangeCellType;   ///< geometry of the smoothing domain

    using TensorD    = Fastor::Tensor< double, nDim >;                   ///< vector of size nDim
    using TensorDD   = Fastor::Tensor< double, nDim, nDim >;             ///< second-order tensor
    using TensorDDD  = Fastor::Tensor< double, nDim, nDim, nDim >;       ///< third-order tensor
    using TensorDDDD = Fastor::Tensor< double, nDim, nDim, nDim, nDim >; ///< fourth-order tensor

    TensorDD _momentsOfInertia_Undeformed; ///< second moments @f$ \boldsymbol{M}^0 @f$ of the undeformed domain
    TensorDD _momentsOfInertia_IntermediateReference; ///< second moments @f$ \boldsymbol{M}^Y @f$ in the
                                                      ///< intermediate reference configuration
    std::array< Eigen::MatrixXd, nDim > _d2N_dYdY;    ///< smoothed second derivatives; `_d2N_dYdY[J](I, B)` is the
                                                      ///< face sum of @f$ \partial N_B/\partial Y_I\,n_J\,dA/V @f$

    /// second gradient of the displacement increment of the last explicit increment, @f$ \partial^2\Delta u_i/
    /// \partial Y_j\partial Y_k @f$; computePhysicsKernelsExplicit() charges its stress gradient to the stabilization
    TensorDDD _d2x_dYdY = TensorDDD( 0.0 );

    /// velocity gradient of the last explicit increment, @f$ \partial\boldsymbol{v}/\partial\boldsymbol{Y} =
    /// \partial\Delta\boldsymbol{u}/\partial\boldsymbol{Y} / \Delta t @f$, for the stabilization of the momentum
    TensorDD _dv_dY = TensorDD( 0.0 );

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
                                                         Fastor::determinant( this->_mp.dY_dX() ) );
    }

    /**
     * @brief The stiffness rows of node @p A: those of GradientEnhancedFiniteStrainParticle plus the stabilization of
     * the displacement rows at the current state, which depends on the second gradient of the displacement increment,
     * @f[
     *   \frac{\partial f^U_{Aj}}{\partial (\partial^2 \Delta u_k / \partial Y_l \partial Y_m)} = D_{im}
     *   \frac{\partial \tau_{ij}}{\partial \Delta F_{kl}}
     * @f]
     * (see momentWeightedSecondGradientOfTestFunction()). The stabilization's dependence on the spatial gradient of
     * the shape functions through the stress gradient of the increment vanishes at a zero increment; it is added in
     * computePhysicsKernels().
     * @param[in] A The node.
     * @param[in] dY_dx The inverse of the gradient of the motion of the particle, @f$ \partial Y / \partial x @f$,
     * the same for all its nodes: computed once per particle by the caller.
     * @param[in,out] rows The derivatives; the stiffness of the node is added.
     */
    void addStiffnessRowsOfNode( int                                         A,
                                 const Fastor::Tensor< double, nDim, nDim >& dY_dx,
                                 typename ParentPointParticle::TangentRows&  rows ) const override;

    /**
     * @brief The smoothed second gradients of the trial shape functions, symmetrized, in the layout of
     * addTangentRowsOfNode() (NSNI::secondGradientsOfNodes()).
     * @return The nDim^2 x nNodes matrix.
     */
    Eigen::MatrixXd secondGradientsOfTrialFunctions() const override
    {
      return NSNI::secondGradientsOfNodes< nDim >( _d2N_dYdY, this->_nNodes );
    }

    /**
     * @brief The second gradient of the displacement increment, @f$ \partial^2\Delta u_i/\partial Y_j\partial Y_k =
     * \Delta q^U_{Bi}\,\partial^2 N_B/\partial Y_j\partial Y_k @f$, with the smoothed second gradients of the shape
     * functions; the part of the interpolation that GradientEnhancedFiniteStrainParticle::interpolateIncrement() does
     * not form.
     * @param[in] dQ Increment of the nodal dofs, node by node.
     * @return @f$ \partial^2\Delta\boldsymbol{u}/\partial\boldsymbol{Y}\partial\boldsymbol{Y} @f$.
     */
    TensorDDD interpolateSecondGradientOfIncrement( const double* dQ ) const;

  public:
    /**
     * @brief Construct the particle and compute the second moments of the undeformed domain.
     *
     * @f$ \boldsymbol{M}^Y @f$ is initialized to @f$ \boldsymbol{M}^0 @f$, since the intermediate reference coincides
     * with the undeformed configuration until the first accepted increment.
     *
     * @param[in] elementID                 Label of the particle.
     * @param[in] nodeCoordinates           Undeformed vertex coordinates, vertex by vertex.
     * @param[in] nNodeCoordiantes          Number of coordinates (unused).
     * @param[in] volume                    Volume (unused, computed from the vertices).
     * @param[in] materialName              Name of a Marmot::Materials::MaterialGradientEnhancedFiniteStrain material.
     * @param[in] materialProperties        Material properties.
     * @param[in] sizeMaterialProperties    Number of material properties.
     * @param[in] approximation             Meshfree approximation used for the shape functions.
     * @param[in] smoothingVolumeUpdateType Update type of the smoothing domain.
     */
    GradientEnhancedFiniteStrainParticleSQCNIxNSNI(
      int                                                     elementID,
      const double*                                           nodeCoordinates,
      int                                                     nNodeCoordiantes,
      double                                                  volume,
      const std::string&                                      materialName,
      const double*                                           materialProperties,
      int                                                     sizeMaterialProperties,
      const Marmot::Discretizations::Meshfree::Approximation& approximation,
      ParentSQCNIParticle::SmoothingDomainUpdateType          smoothingVolumeUpdateType );

    /**
     * @brief Assign the kernel functions and compute the smoothed first and second shape function derivatives.
     *
     * One boundary integration over the current smoothing domain (one point per face) gives both
     * @f$ \partial N_B/\partial\boldsymbol{Y} @f$ (from @f$ N_B @f$ at the face centers) and the second derivatives
     * (from @f$ \partial N_B/\partial\boldsymbol{Y} @f$ at the face centers).
     *
     * @param[in] kernelFunctions Kernel functions of the nodes that support the particle.
     */
    void assignMeshfreeKernelFunctions(
      const std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* >& kernelFunctions ) override
    {
      ParentPointParticle::_assignedKernelFunctions = kernelFunctions;

      ParentPointParticle::_nNodes = static_cast< int >(
        GradientEnhancedFiniteStrainParticle< nDim, nNonlocal >::_assignedKernelFunctions.size() );

      Eigen::Matrix< double, nDim, 1 > coords;
      ParentPointParticle::_mp.getCoordinatesAtCenter( coords.data() );

      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      this->getSmoothingDomainVertexCoordinates( vertexCoordinates.data() );

      ParentPointParticle::_N     = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );
      ParentPointParticle::_dN_dY = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );
      for ( int i = 0; i < nDim; i++ )
        _d2N_dYdY[i] = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         ParentPointParticle::_N.data() );

      Eigen::MatrixXd NBoundary( 1, ParentPointParticle::_nNodes );
      Eigen::MatrixXd dN_dY_Boundary = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

      const LagrangeCellType cell( vertexCoordinates.data(), nDim * nVertices );

      for ( int i = 0; i < cell.getNumberOfFaces(); i++ ) {

        const Eigen::Matrix< double, nDim, 1 > n_dA       = cell.boundarySurfaceVector( i + 1 );
        const Eigen::Matrix< double, nDim, 1 > faceCenter = cell.getFaceCenterCoordinates( i + 1 );

        ParentPointParticle::_meshfreeApproximation
          .computeShapeFunctionsAndGradients( faceCenter.data(),
                                              ParentPointParticle::_assignedKernelFunctions,
                                              NBoundary.data(),
                                              dN_dY_Boundary.data() );

        ParentPointParticle::_dN_dY += n_dA * NBoundary;

        for ( int j = 0; j < nDim; j++ )
          _d2N_dYdY[j] += n_dA( j ) * dN_dY_Boundary;
      }

      const double VSmoothing = cell.volume();
      ParentPointParticle::_dN_dY /= VSmoothing;
      for ( int i = 0; i < nDim; i++ ) {
        _d2N_dYdY[i] /= VSmoothing;
      }

      ParentPointParticle::_T     = ParentPointParticle::_N;
      ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
    }

    /**
     * @brief Residual and tangent of the particle for the increment @p dQ, including the NSNI stabilization.
     *
     * The same as GradientEnhancedFiniteStrainParticle::computePhysicsKernels, plus the stabilization residual and its
     * (approximate) tangent given in the class description. The results are ADDED to @p fInt and @p dFInt_ddQ.
     *
     * @param[in]     dQ        Increment of the nodal dofs since the last accepted state.
     * @param[in,out] fInt      Internal force vector.
     * @param[in,out] dFInt_ddQ Tangent, column-major.
     * @param[in]     timeNew   Time at the end of the increment.
     * @param[in]     dT        Time increment.
     */
    void computePhysicsKernels( const double* dQ, double* fInt, double* dFInt_ddQ, double timeNew, double dT ) override;

    /**
     * @brief Explicit scheme: as GradientEnhancedFiniteStrainParticle::updatePhysicsExplicit, and keep the second
     * gradient and the velocity gradient of the increment for the stabilization (as
     * DisplacementParticleSQCNIxNSNI does).
     * @param[in] dQ      Increment of the nodal dofs, node by node.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT      Time increment.
     */
    void updatePhysicsExplicit( const double* dQ, double timeNew, double dT ) override;

    /**
     * @brief Explicit scheme: the internal force of GradientEnhancedFiniteStrainParticle::computePhysicsKernelsExplicit
     * plus the NSNI stabilization of the displacement rows, @f$ \partial^2 T_A/\partial x_i\partial Y_J\,
     * M^Y_{JK}/J_Y\;\partial\tau_{ij}/\partial Y_K @f$, with the stress gradient of the last increment.
     * @param[in,out] fInt Internal force vector, node by node; the contribution is added.
     */
    void computePhysicsKernelsExplicit( double* fInt ) override;

    /**
     * @brief Lumped momentum of GradientEnhancedFiniteStrainParticle::computeLumpedMomentum plus the NSNI
     * stabilization of the displacement rows, @f$ \rho_0\,\partial\boldsymbol{v}/\partial\boldsymbol{Y}\,
     * \boldsymbol{M}^Y\,\partial T_A/\partial\boldsymbol{Y} / J_Y @f$.
     * @param[in,out] mLumped Lumped momentum vector, node by node; the contribution is added.
     */
    void computeLumpedMomentum( double* mLumped ) const override;

    /**
     * @brief Accept the increment (see GradientEnhancedFiniteStrainParticle::acceptStateAndPosition) and push the
     * second moments forward, @f$ \boldsymbol{M}^Y =
     * J_Y\,\boldsymbol{F}_n\,\boldsymbol{M}^0\,\boldsymbol{F}_n^\mathsf{T}
     * @f$.
     */
    void acceptStateAndPosition() override
    {
      ParentSQCNIParticle::acceptStateAndPosition();

      const auto   F    = this->_mp.dY_dX();
      const double detJ = det( F );

      // I2_ij = FiI * I2_IJ * F_jJ^T * detJ
      _momentsOfInertia_IntermediateReference = F % _momentsOfInertia_Undeformed % transpose( F ) * detJ;
    };
  };

  template < int nDim, int nVertices, int nNonlocal >
  GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::
    GradientEnhancedFiniteStrainParticleSQCNIxNSNI(
      int                                                     elementID,
      const double*                                           vertexCoordinates,
      int                                                     nVertexCoordinates,
      double                                                  volume,
      const std::string&                                      materialName,
      const double*                                           materialProperties,
      int                                                     sizeMaterialProperties,
      const Marmot::Discretizations::Meshfree::Approximation& approximation,
      ParentSQCNIParticle::SmoothingDomainUpdateType          smoothingVolumeUpdateType )
    : GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal >( elementID,
                                                                               vertexCoordinates,
                                                                               nVertexCoordinates,
                                                                               volume,
                                                                               materialName,
                                                                               materialProperties,
                                                                               sizeMaterialProperties,
                                                                               approximation,
                                                                               smoothingVolumeUpdateType )
  {
    // second moments of the undeformed particle domain about its centroid
    // (for an affinely mapped cell this is identical to the former
    // detJ * 16/12 * J J^T quad formula)
    const Eigen::Matrix< double, nDim, nDim > secondMoments = this->_makeUndeformedCell().secondMoments();

    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        _momentsOfInertia_Undeformed( i, j ) = secondMoments( i, j );

    // the intermediate reference coincides with the undeformed one until the first accepted increment;
    // without this, the stabilization of the first increment is scaled by uninitialized memory
    _momentsOfInertia_IntermediateReference = _momentsOfInertia_Undeformed;
  }

  template < int nDim, int nVertices, int nNonlocal >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::computePhysicsKernels(
    const double* dQ,
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

    constexpr int nodeBlockSize = nDim + nNonlocal;

    using NonlocalVector = typename ParentPointParticle::NonlocalVector;

    TensorD                          du;
    TensorDD                         du_dY;
    NonlocalVector                   dn;
    std::array< TensorD, nNonlocal > dn_dY;
    this->interpolateIncrement( dQ, du, du_dY, dn, dn_dY );
    const TensorDDD d2x_dYdY = interpolateSecondGradientOfIncrement( dQ );

    _mp.prepareYourself( timeNew, dT );
    _mp.incrementDeformation( du, du_dY, dn, dn_dY );
    _mp.computeYourself( timeNew, dT );

    const double density0 = _mp.getDensityUndeformed();

    auto v = _mp.getVelocity();
    auto a = _mp.getAcceleration();

    Tensor< double, nDim, nDim > da_ddu( 0.0 );
    Marmot::ContinuumMechanics::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(),
                                                                                 v.data(),
                                                                                 a.data(),
                                                                                 dT,
                                                                                 this->_newmark_beta,
                                                                                 this->_newmark_gamma,
                                                                                 da_ddu.data() );
    _mp.setVelocity( v );
    _mp.setAcceleration( a );

    TensorD r_U( 0.0 );
    double  r_N( 0.0 );

    const auto& S           = _mp.response.S;
    const auto& dLocalField = _mp.response.dL;
    const auto& R           = _mp.response.nonLocalRadius;

    NonlocalVector c;
    for ( int kNl = 0; kNl < nNonlocal; kNl++ )
      c( kNl ) = R( kNl ) * R( kNl );

    const double V0 = this->getVolumeUndeformed();

    const auto& t = _mp.tangents;

    Eigen::Map< Eigen::VectorXd > P( fInt, _nNodes * nodeBlockSize );
    Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

    const auto dY_dx = evaluate( inv( _mp.dx_dY() ) );

    // Gradient of the Kirchhoff stress along Y, which is what the NSNI stabilization charges.
    // Only the deformation-gradient term survives here: the micropolar reference also carried
    // dS_dW and dS_ddWdY, and like it we do NOT include the dS_dN contribution of the nonlocal
    // field.  Note that the tangent of the stabilization below is APPROXIMATE: it differentiates
    // d2x_dYdY but not dS_dDeltaF itself, i.e. it omits d2tau/dF2 (not exposed by the material
    // interface) and the dependence of dS_dDeltaF on the nonlocal field (the damage). The nonlocal
    // rows are exact; see the module test for the measured size of the omitted terms.
    const auto dS_dY = evaluate( einsum< ijmM, mMK >( t.dS_dDeltaF, d2x_dYdY ) );

    const Eigen::MatrixXd d2N_dYdY = secondGradientsOfTrialFunctions();

    typename ParentPointParticle::TangentRows rows;
    for ( int A = 0; A < _nNodes; A++ ) {

      const double  T_A     = _T( A );
      const auto    dT_A_dY = TensorMap< const double, nDim >( _dT_dY.col( A ).data() );
      const TensorD dT_A_dx = einsum< ji, j >( dY_dx, dT_A_dY );
      const TensorD dT_A_dX = einsum< ji, j >( _mp.dY_dX(), dT_A_dY );

      const int idxA_u = nodeBlockSize * A;
      const int idxA_n = nodeBlockSize * A + nDim;

      r_U = ( +einsum< i, ij >( dT_A_dx, S ) ) * V0;

      // add inertia
      r_U += density0 * a * T_A * V0;

      const TensorDD d2NA_dxdY_x_MOIScaled = momentWeightedSecondGradientOfTestFunction( A, dY_dx );

      TensorD rU_Stab = einsum< iK, ijK >( d2NA_dxdY_x_MOIScaled, dS_dY );

      r_U += rU_Stab;

      P.template segment< nDim >( idxA_u ) += Eigen::Map< Eigen::Matrix< double, nDim, 1 > >( r_U.data() );

      for ( int kNl = 0; kNl < nNonlocal; kNl++ ) {
        const TensorD dn_dX = einsum< ji, j >( _mp.dY_dX(), dn_dY[kNl] );

        r_N = evaluate( ( T_A * dn( kNl ) + c( kNl ) * einsum< i, i >( dT_A_dX, dn_dX ) - T_A * dLocalField( kNl ) ) *
                        V0 )
                .toscalar();

        P( idxA_n + kNl ) += r_N;
      }

      rows.setZero();
      addStiffnessRowsOfNode( A, dY_dx, rows );
      rows.df_dIncrement
        .template topLeftCorner< nDim,
                                 nDim >() += density0 * T_A * V0 *
                                             Eigen::Map< const Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >(
                                               da_ddu.data() );

      // the stabilization through the stress gradient of the increment: D_kK dTau_ij/dY_K dN_B/dx_i
      NSNI::addStressGradientRows< nDim >( d2NA_dxdY_x_MOIScaled, dS_dY, dY_dx, rows );

      addTangentRowsOfNode( rows, A, _N, _dN_dY, d2N_dYdY, K );
    }
  }

  template < int nDim, int nVertices, int nNonlocal >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::addStiffnessRowsOfNode(
    int                                         A,
    const Fastor::Tensor< double, nDim, nDim >& dY_dx,
    typename ParentPointParticle::TangentRows&  rows ) const
  {
    ParentPointParticle::addStiffnessRowsOfNode( A, dY_dx, rows );
    NSNI::addSecondGradientRows< nDim >( momentWeightedSecondGradientOfTestFunction( A, dY_dx ),
                                         this->_mp.tangents.dS_dDeltaF,
                                         rows );
  }

  template < int nDim, int nVertices, int nNonlocal >
  typename GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::TensorDDD GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
    nDim,
    nVertices,
    nNonlocal >::interpolateSecondGradientOfIncrement( const double* dQ ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    constexpr int nodeBlockSize = nDim + nNonlocal;

    TensorDDD d2x_dYdY( 0.0 );
    for ( int B = 0; B < this->_nNodes; B++ )
      d2x_dYdY += einsum< i, jk >( TensorD( dQ + nodeBlockSize * B ),
                                   NSNI::secondGradientOfNode< nDim >( _d2N_dYdY, B ) );
    return d2x_dYdY;
  }

  template < int nDim, int nVertices, int nNonlocal >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::updatePhysicsExplicit(
    const double* dQ,
    double        timeNew,
    double        dT )
  {
    _d2x_dYdY = interpolateSecondGradientOfIncrement( dQ );

    ParentPointParticle::updatePhysicsExplicit( dQ, timeNew, dT );

    if ( dT <= 0.0 )
      return;

    TensorDD I;
    I.eye();
    _dv_dY = ( this->_mp.dx_dY() - I ) / dT;
  }

  template < int nDim, int nVertices, int nNonlocal >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::computePhysicsKernelsExplicit(
    double* fInt )
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;
    using ijmM = Index< i_, j_, m_, M_ >;
    using mMK  = Index< m_, M_, K_ >;
    using ijK  = Fastor::Index< i_, j_, K_ >;

    ParentPointParticle::computePhysicsKernelsExplicit( fInt );

    constexpr int nodeBlockSize = nDim + nNonlocal;

    const auto     dS_dY = evaluate( einsum< ijmM, mMK >( this->_mp.tangents.dS_dDeltaF, _d2x_dYdY ) );
    const TensorDD dY_dx = inv( this->_mp.dx_dY() );

    Eigen::Map< Eigen::VectorXd > P( fInt, this->_nNodes * nodeBlockSize );

    for ( int A = 0; A < this->_nNodes; A++ ) {
      const TensorDD d2NA_dxdY_x_MOIScaled = momentWeightedSecondGradientOfTestFunction( A, dY_dx );
      const TensorD  rU_Stab               = einsum< iK, ijK >( d2NA_dxdY_x_MOIScaled, dS_dY );

      P.template segment< nDim >( nodeBlockSize *
                                  A ) += Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( rU_Stab.data() );
    }
  }

  template < int nDim, int nVertices, int nNonlocal >
  void GradientEnhancedFiniteStrainParticleSQCNIxNSNI< nDim, nVertices, nNonlocal >::computeLumpedMomentum(
    double* mLumped ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    ParentPointParticle::computeLumpedMomentum( mLumped );

    constexpr int nodeBlockSize = nDim + nNonlocal;

    const double   density0         = this->_mp.getDensityUndeformed();
    const double   detJIntermediate = determinant( this->_mp.dY_dX() );
    const TensorDD dv_dY_x_Y2       = einsum< ij, jk >( _dv_dY, _momentsOfInertia_IntermediateReference );

    for ( int A = 0; A < this->_nNodes; A++ ) {
      const TensorD dT_A_dY = TensorD( this->_dT_dY.col( A ).data() );
      const TensorD aux     = einsum< ij, j >( dv_dY_x_Y2, dT_A_dY );
      for ( int i = 0; i < nDim; i++ )
        mLumped[nodeBlockSize * A + i] += density0 * aux( i ) / detJIntermediate;
    }
  }

} // namespace Marmot::Discretizations::Meshfree
