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

#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeQuadHexCell.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <cmath>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::GradientEnhancedFiniteStrainParticleSQCNI
   * @brief Stabilized quasi-conforming nodal integration (SQCNI / SNNI) particle for gradient-enhanced finite-strain
   *        materials without a micropolar continuum.
   *
   * The particle is a quadrilateral (2D) or hexahedral (3D) smoothing domain with one material point at its centroid.
   * The weak forms, the material point and the tangent are those of GradientEnhancedFiniteStrainParticle; what
   * changes is how the shape function gradients are obtained:
   *
   * - **Smoothed gradients.** As in stabilized conforming nodal integration (Chen et al., 2001), the gradient is
   *   replaced by its average over the smoothing domain @f$ \Omega_Y @f$, turned into a boundary integral by the
   *   divergence theorem and evaluated with one point per face @f$ f @f$ (the face center @f$ \boldsymbol{Y}_f @f$),
   *   @f[
   *     \frac{\partial N_B}{\partial Y_i} \approx \frac{1}{V_{\Omega_Y}} \sum_f N_B(\boldsymbol{Y}_f)\,n_i\,dA_f .
   *   @f]
   *   The shape function values @f$ N_B @f$ remain the point values at the center. Both the displacement gradient
   *   @f$ \Delta\boldsymbol{F} @f$ and the gradient of the nonlocal field use the smoothed gradients.
   * - **Faces.** The particle has vertices and faces and therefore supports the distributed loads `PRESSURE` and
   *   `CWFCORRECTION` (computeDistributedLoad()); the point particle has none.
   *
   * The smoothing domain is stored by the displacements of its vertices (appended to the state variables as
   * `vertex displacements`) and is updated at each accepted increment from the deformation gradient
   * @f$ \boldsymbol{F}_n @f$ of the material point, according to SmoothingDomainUpdateType:
   * `DeformationGradient` (SQCNI: the domain follows @f$ \boldsymbol{F}_n @f$), `None` (SNNI: the domain is only
   * translated with the particle), `RotationOnly` (the rotation @f$ \boldsymbol{R} @f$ of the polar decomposition)
   * and `RotationAndPrincipalStretch` (@f$ \boldsymbol{R} @f$ times the diagonal of
   * @f$ \boldsymbol{R}^\mathsf{T}\boldsymbol{F}_n @f$). Only for `DeformationGradient` does the smoothing domain
   * coincide with the physical domain of the particle.
   *
   * @tparam nDim      Spatial dimension (2: plane strain, 3: 3D).
   * @tparam nVertices Number of vertices of the smoothing domain (4: Quad, 8: Hexa).
   */
  template < int nDim, int nVertices >
  class GradientEnhancedFiniteStrainParticleSQCNI : public GradientEnhancedFiniteStrainParticle< nDim > {

    using TensorD          = Fastor::Tensor< double, nDim >;   ///< vector of size nDim
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >; ///< coordinate vector

  public:
    /// @brief How the smoothing domain follows the deformation at each accepted increment.
    enum SmoothingDomainUpdateType {
      None,                       ///< translation only (SNNI)
      DeformationGradient,        ///< mapped by @f$ \boldsymbol{F}_n @f$ (SQCNI)
      RotationOnly,               ///< mapped by the rotation @f$ \boldsymbol{R} @f$ of @f$ \boldsymbol{F}_n @f$
      RotationAndPrincipalStretch ///< mapped by @f$
                                  ///< \boldsymbol{R}\,\mathrm{diag}(\boldsymbol{R}^\mathsf{T}\boldsymbol{F}_n) @f$
    };

  protected:
    using LagrangeCellType = MarmotLagrangeCell< nDim, nVertices >; ///< geometry of the smoothing domain

    const SmoothingDomainUpdateType _smoothingVolumeUpdateType;     ///< update type of the smoothing domain

    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed; ///< undeformed vertex coordinates

    double*
      _vertexDisplacements_SmoothingDomain; ///< vertex displacements of the smoothing domain (in the state vector)

    using ParentPointParticle = GradientEnhancedFiniteStrainParticle< nDim >; ///< the point particle base

    /// \brief Build a Lagrange cell from the current (deformed) smoothing domain vertex coordinates
    /// \return The smoothing domain in the intermediate reference configuration.
    LagrangeCellType _makeSmoothingDomainCell() const
    {
      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getVertexCoordinates( vertexCoordinates.data() );
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices );
    }

    /// \brief Build a Lagrange cell from the undeformed vertex coordinates
    /// \return The smoothing domain in the undeformed configuration.
    LagrangeCellType _makeUndeformedCell() const
    {
      return LagrangeCellType( _vertexCoordinates_Undeformed.data(), nDim * nVertices );
    }

  public:
    /**
     * @brief Vertex coordinates of the smoothing domain: undeformed coordinates plus `vertex displacements`.
     * @param[out] coordinates Coordinates, vertex by vertex (nDim x nVertices).
     */
    virtual void getVertexCoordinates( double* coordinates ) const override;

    /**
     * @brief Vertex coordinates for visualization, identical to getVertexCoordinates().
     * @param[out] coordinates Coordinates, vertex by vertex.
     */
    virtual void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getVertexCoordinates( coordinates );
    };

    /**
     * @brief Number of vertices.
     * @return nVertices.
     */
    virtual int getNumberOfVertices() const override { return nVertices; };

    /**
     * @brief Shape of the particle.
     * @return Shape name of the undeformed smoothing domain cell.
     */
    virtual std::string getParticleShape() const override { return _makeUndeformedCell().getCellShape(); }

    /**
     * @brief Construct the particle from its vertices.
     *
     * The material point is placed at the centroid of the vertices, with the volume of the cell they span; the
     * @p volume argument is not used.
     *
     * @param[in] elementID                 Label of the particle.
     * @param[in] nodeCoordinates           Undeformed vertex coordinates, vertex by vertex.
     * @param[in] nNodeCoordiantes          Number of coordinates (unused).
     * @param[in] volume                    Volume (unused, computed from the vertices).
     * @param[in] materialName              Name of a MarmotMaterialGradientEnhancedFiniteStrain material.
     * @param[in] materialProperties        Material properties.
     * @param[in] sizeMaterialProperties    Number of material properties.
     * @param[in] approximation             Meshfree approximation used for the shape functions.
     * @param[in] smoothingVolumeUpdateType Update type of the smoothing domain.
     */
    GradientEnhancedFiniteStrainParticleSQCNI( int                                elementID,
                                               const double*                      nodeCoordinates,
                                               int                                nNodeCoordiantes,
                                               double                             volume,
                                               const std::string&                 materialName,
                                               const double*                      materialProperties,
                                               int                                sizeMaterialProperties,
                                               const MarmotMeshfreeApproximation& approximation,
                                               const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    /**
     * @brief Assign the kernel functions and compute the smoothed shape function gradients.
     *
     * @f$ N_B @f$ is evaluated at the material point position, @f$ \partial N_B/\partial\boldsymbol{Y} @f$ as the
     * boundary integral over the current smoothing domain with one point per face (see the class description); the
     * test functions are set equal to the trial functions (before any VCI correction).
     *
     * @param[in] kernelFunctions Kernel functions of the nodes that support the particle.
     */
    virtual void assignMeshfreeKernelFunctions(
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) override;

    /// \brief Get the smoothing volume of the particle
    /// \return The smoothing volume of the particle
    /// \details The smoothing volume is computed as the volume of the element in an updated configuration
    ///         that is obtained by applying (parts of) the deformation gradient to the element in the undeformed
    ///         configuration. In general, this is not consistent with the actual deformation of the particle.
    ///
    virtual double getSmoothingVolume() const { return _makeSmoothingDomainCell().volume(); }

    /**
     * @brief Center coordinates: the material point position of the last accepted state.
     * @param[out] coordinates Coordinates (nDim values).
     */
    virtual void getCenterCoordinates( double* coordinates ) const override
    {
      this->_mp.getCoordinatesAtCenter( coordinates );
    }

    /**
     * @brief Centroid of a cell spanned by given vertices.
     * @param[in] vertexCoordinates Vertex coordinates.
     * @return The centroid.
     */
    virtual CoordinatesSized getCenterFromVertices(
      const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices ).centroid();
    }

    /**
     * @brief Volume of a cell spanned by given vertices.
     * @param[in] vertexCoordinates Vertex coordinates.
     * @return The volume.
     */
    virtual double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices ).volume();
    }

    /**
     * @brief Supported distributed loads, on the faces of the particle domain.
     * @return `PRESSURE` and `CWFCORRECTION`.
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      using Parent = GradientEnhancedFiniteStrainParticle< nDim >;
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes =
        { { "PRESSURE", Parent::Pressure }, { "CWFCORRECTION", Parent::CWFCorrection } };
      return _supportedDistributedLoadTypes;
    };

    /**
     * @brief Volume of the last accepted state, @f$ V_0\det\boldsymbol{F}_n @f$.
     * @return The volume in the intermediate reference configuration (the current increment is not included).
     */
    virtual double getVolumeDeformed() const
    {
      return this->getVolumeUndeformed() * Fastor::determinant( this->_mp.dY_dX() );
    }

    /**
     * @brief Number of state variables: those of the point particle plus the vertex displacements.
     * @return Required size of the state variable vector.
     */
    virtual int getNumberOfRequiredStateVars() const override
    {
      return GradientEnhancedFiniteStrainParticle< nDim >::getNumberOfRequiredStateVars() + nDim * nVertices;
    };

    /**
     * @brief Assign the state variable vector; the last nDim x nVertices entries hold the vertex displacements.
     * @param[in,out] stateVars  State variable vector.
     * @param[in]     nStateVars Its size.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      GradientEnhancedFiniteStrainParticle< nDim >::assignStateVars( stateVars, nStateVars - nDim * nVertices );
      _vertexDisplacements_SmoothingDomain = stateVars + nStateVars - nDim * nVertices;
    }

    /**
     * @brief Access a state variable; `vertex displacements` are those of the smoothing domain.
     * @param[in] stateName Name of the state variable.
     * @param[in] qp        Evaluation point (unused).
     * @return View on the state variable.
     */
    virtual StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" ) {
        return { _vertexDisplacements_SmoothingDomain, nDim * nVertices };
      }
      return GradientEnhancedFiniteStrainParticle< nDim >::getStateView( stateName, qp );
    }

    /**
     * @brief Distributed loads on a face of the particle.
     *
     * Both loads use the surface vector @f$ \boldsymbol{N}\,dA_Y @f$ and center @f$ \boldsymbol{Y}_f @f$ of the face
     * in the intermediate reference configuration (getIntermediateConfBoundaryVector()), test functions
     * @f$ T_A(\boldsymbol{Y}_f) @f$ evaluated by the meshfree approximation at the face center, and are subtracted
     * from @p fExt (sign convention of the internal force), with the consistent tangent with respect to the
     * displacement increment:
     *
     * - `PRESSURE` (follower pressure @f$ p @f$ = `load[0]`):
     *   @f$ \boldsymbol{f} = \Delta J\,\Delta\boldsymbol{F}^{-\mathsf{T}}\,p\,\boldsymbol{N}\,dA_Y @f$,
     *   @f$ f_{Ai} \mathrel{-}= T_A(\boldsymbol{Y}_f)\,f_i @f$.
     * - `CWFCORRECTION` (consistent weak form correction): the boundary term of the momentum weak
     *   form, @f$ \boldsymbol{t} = \boldsymbol{\tau}\,\Delta\boldsymbol{F}^{-\mathsf{T}}\,\boldsymbol{N}\,dA_Y / J_Y
     *   @f$ with @f$ J_Y = \det\boldsymbol{F}_n @f$, @f$ f_{Ai} \mathrel{-}= T_A(\boldsymbol{Y}_f)\,t_i @f$; its
     * tangent contains the geometric part from @f$ \Delta\boldsymbol{F}^{-\mathsf{T}} @f$, the material part from
     *   @f$ \partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F} @f$ and the coupling
     *   @f$ \partial\boldsymbol{\tau}/\partial\bar{N} @f$ to the nonlocal dofs. `load[0]` (or no load) selects
     *   the corrected components as a bit mask: 1 = x, 2 = y, 4 = z, 0 = all.
     *
     * Neither load contributes to the nonlocal residual. Both read the state of the last computePhysicsKernels().
     *
     * @param[in]     type      Load type (see getSupportedDistributedLoadTypes()).
     * @param[in]     surfaceID Face id of the smoothing domain (1-based).
     * @param[in]     load      Load values (`PRESSURE`: the pressure; `CWFCORRECTION`: component mask, may be null).
     * @param[in,out] fExt      Load vector.
     * @param[in,out] dExt_dQ   Load tangent, column-major.
     * @param[in]     timeNew   Time at the end of the increment (unused).
     * @param[in]     dT        Time increment (unused).
     * @throws std::invalid_argument for an unknown load type.
     */
    virtual void computeDistributedLoad( int           type,
                                         int           surfaceID,
                                         const double* load,
                                         double*       fExt,
                                         double*       dExt_dQ,
                                         double        timeNew,
                                         double        dT ) const override;

    /**
     * @brief Surface vector and center of a face in the intermediate reference configuration.
     *
     * For `DeformationGradient` they are taken from the current smoothing domain, which coincides with the
     * particle. For the other update types the undeformed surface vector is mapped by Nanson's formula,
     * @f$ \boldsymbol{N}\,dA_Y = J_n\,\boldsymbol{F}_n^{-\mathsf{T}}\,\boldsymbol{N}\,dA_0 @f$, and the undeformed
     * face center is translated by the particle displacement.
     *
     * @param[in] boundaryFaceID Face id (1-based).
     * @return Tuple of the surface vector @f$ \boldsymbol{N}\,dA_Y @f$ and the face center @f$ \boldsymbol{Y}_f @f$.
     */
    std::tuple< TensorD, TensorD > getIntermediateConfBoundaryVector( int boundaryFaceID ) const;

    /**
     * @brief Add the boundary term @f$ T_A(\boldsymbol{Y}_f)\,P_C(\boldsymbol{Y}_f)\,N_i\,dA_Y @f$ of the VCI
     * integration constraint for a face of the particle.
     *
     * Overrides the point version: the surface vector and the location come from getIntermediateConfBoundaryVector(),
     * and @f$ T_A @f$ and @f$ P_C @f$ are evaluated at the face center.
     *
     * @param[in,out] R_AiC_RowMajor        Constraint residual @f$ R_{AiC} @f$, row-major.
     * @param[in]     boundarySurfaceVector Unused.
     * @param[in]     boundaryFaceID        Face id (1-based).
     */
    virtual void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                                      const double* boundarySurfaceVector,
                                                      int           boundaryFaceID )
    {
      using namespace Fastor;

      const auto [N_dAY, Y_N]   = getIntermediateConfBoundaryVector( boundaryFaceID );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         TBoundary.data() );

      // get P for the exact integration location at the boundary
      auto PBoundary = ParentPointParticle::_P;
      Math::computeMonomialBasis( ParentPointParticle::_vciOrder, CoordinatesSized( Y_N.data() ), PBoundary );

      for ( int A = 0; A < ParentPointParticle::_nNodes; A++ )
        for ( int i = 0; i < nDim; i++ )
          for ( int C = 0; C < ParentPointParticle::_nVCIConstraints; C++ )
            R_AiC_RowMajor[A * ( nDim * ParentPointParticle::_nVCIConstraints ) +
                           i * ParentPointParticle::_nVCIConstraints + C] += TBoundary( A ) * PBoundary( C ) * N_dAY[i];
    };

    /**
     * @brief Coordinates of the evaluation points: the face centers of the smoothing domain.
     * @param[out] coordinates Coordinates, face by face (nDim x number of faces).
     */
    virtual void getEvaluationCoordinates( double* coordinates ) const
    {
      const auto cell = _makeSmoothingDomainCell();

      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > faceCenters( coordinates,
                                                                               nDim,
                                                                               cell.getNumberOfFaces() );

      for ( int i = 0; i < cell.getNumberOfFaces(); i++ )
        faceCenters.col( i ) = cell.getFaceCenterCoordinates( i + 1 );
    }

    /**
     * @brief Center of a face of the current smoothing domain.
     * @param[in]  faceID      Face id (1-based).
     * @param[out] coordinates Face center (nDim values).
     */
    virtual void getFaceCoordinates( int faceID, double* coordinates ) const
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > faceCenter( coordinates );
      faceCenter = _makeSmoothingDomainCell().getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Number of evaluation points.
     * @return Number of faces of the smoothing domain.
     */
    virtual int getNumberOfEvaluationPoints() const { return _makeUndeformedCell().getNumberOfFaces(); };

  private:
    /// @brief Update the vertex displacements of the smoothing domain from @f$ \boldsymbol{F}_n @f$ and the particle
    /// displacement, according to the SmoothingDomainUpdateType.
    void _updateVertexDisplacementsFromMaterialPointDeformation();

    /// @brief Update the smoothing domain, then set the center to the accepted material point position.
    virtual void updateParticlePositionToReferenceIntermediate() override
    {
      _updateVertexDisplacementsFromMaterialPointDeformation();

      getCenterCoordinates( ParentPointParticle::_centerReferenceIntermediate.data() );
    };
  };

  template < int nDim, int nVertices >
  GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::GradientEnhancedFiniteStrainParticleSQCNI(
    int                                                  elementID,
    const double*                                        vertexCoordinates,
    int                                                  nVertexCoordinates,
    double                                               volume,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation,
    const SmoothingDomainUpdateType                      smoothingVolumeUpdateType )
    : GradientEnhancedFiniteStrainParticle< nDim >( elementID,
                                                    getCenterFromVertices(
                                                      Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                        vertexCoordinates ) )
                                                      .data(),
                                                    CoordinatesSized::RowsAtCompileTime,
                                                    getVolumeFromVertices(
                                                      Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >(
                                                        vertexCoordinates ) ),
                                                    materialName,
                                                    materialProperties,
                                                    sizeMaterialProperties,
                                                    approximation ),
      _smoothingVolumeUpdateType( smoothingVolumeUpdateType ),
      _vertexCoordinates_Undeformed( vertexCoordinates )

  {
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::computeDistributedLoad( int type,
                                                                                             int boundaryFaceID,
                                                                                             const double* load_,
                                                                                             double*       fExt,
                                                                                             double*       dFExt_ddQ,
                                                                                             double        timeNew,
                                                                                             double        dT ) const
  {

    switch ( type ) {

    case GradientEnhancedFiniteStrainParticle< nDim >::Pressure: {

      const auto&   _mp           = GradientEnhancedFiniteStrainParticle< nDim >::_mp;
      const auto&   _nNodes       = GradientEnhancedFiniteStrainParticle< nDim >::_nNodes;
      constexpr int nodeBlockSize = nDim + 1;

      const auto [N_dAY, Y_N] = getIntermediateConfBoundaryVector( boundaryFaceID );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye;
      Eye.eye();

      // apply Nanson's formula
      const auto deltaF = _mp.dx_dY();

      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;

      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );

      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );

      Eigen::MatrixXd testBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        /* const double T_A = GradientEnhancedFiniteStrainParticle< nDim >::_T( A ); */
        r_U = testBoundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = TensorMap< const double, nDim >(
            GradientEnhancedFiniteStrainParticle< nDim >::_dN_dY.col( B ).data() );

          const Tensor< double, nDim, nDim > df_ddQU_B = testBoundary( A ) * einsum< ijk, k >( df_dDeltaF, dN_B_dY );

          {
            using namespace Eigen;
            K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< Matrix< double, nDim, nDim > >(
              torowmajor( df_ddQU_B ).data() );
          }
        }
      }

      break;
    }
    case GradientEnhancedFiniteStrainParticle< nDim >::CWFCorrection: {
      const auto&   _mp           = GradientEnhancedFiniteStrainParticle< nDim >::_mp;
      const auto&   _nNodes       = GradientEnhancedFiniteStrainParticle< nDim >::_nNodes;
      constexpr int nodeBlockSize = nDim + 1;

      const auto [N_dAY, Y_N] = getIntermediateConfBoundaryVector( boundaryFaceID );

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;

      const auto& S = _mp.response.S;
      const auto& t = _mp.tangents;

      // the internal force integrates tau over the undeformed volume V0 = V_Y / J_Y, so the boundary term of the same
      // weak form is tau * deltaF^-T * N dA_Y / J_Y, with J_Y the Jacobian of the intermediate configuration
      const auto                         deltaF    = _mp.dx_dY();
      const Tensor< double, nDim, nDim > deltaFInv = inverse( deltaF );
      const double                       JY        = determinant( _mp.dY_dX() );

      const TensorD v        = transpose( deltaFInv ) % N_dAY / JY;
      TensorD       traction = S % v;

      // d traction_i / d deltaF_mM: the geometric part from d(deltaF^-T), the material part from d tau
      Tensor< double, nDim, nDim, nDim > dTraction_dDeltaF( 0.0 );
      for ( int i = 0; i < nDim; ++i )
        for ( int m = 0; m < nDim; ++m )
          for ( int M = 0; M < nDim; ++M )
            for ( int j = 0; j < nDim; ++j )
              dTraction_dDeltaF( i, m, M ) += -S( i, j ) * v( m ) * deltaFInv( M, j ) +
                                              t.dS_dDeltaF( i, j, m, M ) * v( j );

      TensorD dTraction_dN = t.dS_dN % v;

      // load[0] selects the corrected components as a bit mask (1 = x, 2 = y, 4 = z; 0 or no load = all): on a face
      // where only the normal displacement is prescribed (symmetry plane, frictionless platen) only that component of
      // the traction is a reaction, the tangential one is a natural (zero) condition
      const int mask = load_ ? static_cast< int >( std::lround( load_[0] ) ) : 0;
      TensorD   sel;
      for ( int i = 0; i < nDim; ++i )
        sel( i ) = ( mask == 0 || ( mask >> i ) & 1 ) ? 1.0 : 0.0;
      for ( int i = 0; i < nDim; ++i ) {
        traction( i ) *= sel( i );
        dTraction_dN( i ) *= sel( i );
        for ( int m = 0; m < nDim; ++m )
          for ( int M = 0; M < nDim; ++M )
            dTraction_dDeltaF( i, m, M ) *= sel( i );
      }

      Eigen::RowVectorXd testBoundary = Eigen::RowVectorXd::Zero( ParentPointParticle::_nNodes );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( Y_N.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         testBoundary.data() );

      for ( int A = 0; A < _nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;

        const TensorD r_U = testBoundary( A ) * traction;
        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< const Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < _nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const int  idxB_n  = nodeBlockSize * B + nDim;
          const auto dN_B_dY = TensorMap< const double, nDim >(
            GradientEnhancedFiniteStrainParticle< nDim >::_dN_dY.col( B ).data() );
          const double N_B = GradientEnhancedFiniteStrainParticle< nDim >::_N( B );

          Tensor< double, nDim, nDim > k_UU( 0.0 );
          for ( int i = 0; i < nDim; ++i )
            for ( int m = 0; m < nDim; ++m )
              for ( int M = 0; M < nDim; ++M )
                k_UU( i, m ) += testBoundary( A ) * dTraction_dDeltaF( i, m, M ) * dN_B_dY( M );

          const TensorD k_UN = testBoundary( A ) * N_B * dTraction_dN;

          {
            using namespace Eigen;
            K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< const Matrix< double, nDim, nDim > >(
              torowmajor( k_UU ).data() );
            K.template block< nDim, 1 >( idxA_u, idxB_n ) -= Map< const Matrix< double, nDim, 1 > >( k_UN.data() );
          }
        }
      }
      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices >
  std::tuple< typename GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::TensorD,
              typename GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::TensorD >
  GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::getIntermediateConfBoundaryVector(
    int boundaryFaceID ) const
  {

    TensorD N_dAY;
    TensorD Y;

    if ( _smoothingVolumeUpdateType == DeformationGradient ) {
      // For the full SQCNI case, we actually operate on the deformed element.
      // This means, that the deformed smoothing domain is consistent with the physical domain of the particle.
      const auto cell = _makeSmoothingDomainCell();

      const Eigen::Matrix< double, nDim, 1 > n_dA = cell.boundarySurfaceVector( boundaryFaceID );
      const Eigen::Matrix< double, nDim, 1 > y    = cell.getFaceCenterCoordinates( boundaryFaceID );

      for ( int i = 0; i < nDim; i++ ) {
        N_dAY[i] = n_dA[i];
        Y[i]     = y[i];
      }
    }
    else {
      // In all other cases, the deformed smoothing domain is not consistent with the physical domain of the particle.
      // Accordingly, it is not reasonable to let the boundary vector be acting on the deformed smoothing domain.
      // Instead, we compute the boundary vector in the undeformed configuration and apply the deformation gradient
      // using Nanson's formula. Also, the let the origin of the boundary segment be the center of the particle in the
      // deformed configuration.

      const auto cellUndeformed = _makeUndeformedCell();

      const Eigen::Matrix< double, nDim, 1 > N_dA0_eigen = cellUndeformed.boundarySurfaceVector( boundaryFaceID );
      const Eigen::Matrix< double, nDim, 1 > Y0          = cellUndeformed.getFaceCenterCoordinates( boundaryFaceID );

      const auto FY = this->_mp.dY_dX();

      // apply Nanson's formula to shift to Y configuration:
      using namespace Fastor;
      N_dAY = determinant( FY ) * transpose( inverse( FY ) ) % TensorD( N_dA0_eigen.data() );

      // and the center of the boundary segment in the deformed configuration
      this->_mp.getCenterDisplacement( Y.data() );
      Y += TensorD( Y0.data() );
    }

    return { N_dAY, Y };
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::getVertexCoordinates( double* coordinates ) const
  {

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesDeformed( coordinates );
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacementsSmoothingDomain(
      _vertexDisplacements_SmoothingDomain );

    coordinatesDeformed = _vertexCoordinates_Undeformed + vertexDisplacementsSmoothingDomain;
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim,
                                                  nVertices >::_updateVertexDisplacementsFromMaterialPointDeformation()
  {

    Eigen::Matrix< double, nDim, nDim > F;

    switch ( _smoothingVolumeUpdateType ) {

    case SmoothingDomainUpdateType::None: F.setIdentity(); break;

    case SmoothingDomainUpdateType::DeformationGradient: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );
      break;
    }

    case SmoothingDomainUpdateType::RotationOnly: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      F = R;
      break;
    }

    case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim >::_mp.dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      Eigen::JacobiSVD< Eigen::MatrixXd > svd;
      svd.compute( F, Eigen::ComputeFullU | Eigen::ComputeFullV );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      Eigen::Matrix< double, nDim, nDim > U_ = Eigen::Matrix< double, nDim, nDim >::Identity();
      U_.diagonal()                          = ( R.transpose() * F ).diagonal();

      F = R * U_;

      break;
    }
    }

    Eigen::Matrix< double, nDim, nVertices > coordinatesRelToCenter0 = _vertexCoordinates_Undeformed;
    for ( int i = 0; i < nDim; i++ ) {
      coordinatesRelToCenter0.row( i ).array() -= this->_centerCoordinatesUndeformed( i );
    }

    Eigen::Matrix< double, nDim, nVertices > vertexCoordinatesDeformed = F * coordinatesRelToCenter0;

    Eigen::Matrix< double, nDim, 1 > _mp_displacement;
    GradientEnhancedFiniteStrainParticle< nDim >::_mp.getCenterDisplacement( _mp_displacement.data() );

    for ( int i = 0; i < nDim; i++ ) {
      vertexCoordinatesDeformed.row( i ).array() += this->_centerCoordinatesUndeformed( i ) + _mp_displacement( i );
    }

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacements( _vertexDisplacements_SmoothingDomain );
    vertexDisplacements = vertexCoordinatesDeformed - _vertexCoordinates_Undeformed;
  }

  template < int nDim, int nVertices >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices >::assignMeshfreeKernelFunctions(
    const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions )
  {

    ParentPointParticle::_assignedKernelFunctions = kernelFunctions;

    ParentPointParticle::_nNodes = GradientEnhancedFiniteStrainParticle< nDim >::_assignedKernelFunctions.size();

    Eigen::Matrix< double, nDim, 1 > coords;
    ParentPointParticle::_mp.getCoordinatesAtCenter( coords.data() );

    ParentPointParticle::_N     = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );
    ParentPointParticle::_dN_dY = Eigen::MatrixXd::Zero( nDim, ParentPointParticle::_nNodes );

    ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                       ParentPointParticle::_assignedKernelFunctions,
                                                                       ParentPointParticle::_N.data() );

    Eigen::MatrixXd NBoundary( 1, ParentPointParticle::_nNodes );

    const auto cell = _makeSmoothingDomainCell();

    for ( int i = 0; i < cell.getNumberOfFaces(); i++ ) {

      const Eigen::Matrix< double, nDim, 1 > n_dA       = cell.boundarySurfaceVector( i + 1 );
      const Eigen::Matrix< double, nDim, 1 > faceCenter = cell.getFaceCenterCoordinates( i + 1 );

      ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( faceCenter.data(),
                                                                         ParentPointParticle::_assignedKernelFunctions,
                                                                         NBoundary.data() );

      ParentPointParticle::_dN_dY += n_dA * NBoundary;
    }

    ParentPointParticle::_dN_dY /= cell.volume();

    ParentPointParticle::_T     = ParentPointParticle::_N;
    ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
  }

} // namespace Marmot::Meshfree
