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
#include "Marmot/MarmotParticleSectionIntegrals.h"
#include <Eigen/Core>
#include <Eigen/Dense>
/* Fastor/Fastor.h ends with `#pragma warning( default : ... )`; scope it so it neither overrides the build's warning
 * flags nor leaks into consumers. */
#ifdef _MSC_VER
#  pragma warning( push )
#endif
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#  pragma warning( pop )
#endif
#include <array>
#include <cmath>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::GradientEnhancedFiniteStrainParticleSQCNI
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
   * The particle geometry and the smoothing domain are stored by the displacements of their vertices (appended to the
   * state variables as `vertex displacements` and `smoothing vertex displacements`). The geometry follows
   * @f$ \boldsymbol{F}_n @f$; its vertices are getVertexCoordinates(). The smoothing domain is updated at each accepted
   * increment from the deformation gradient
   * @f$ \boldsymbol{F}_n @f$ of the material point, according to SmoothingDomainUpdateType:
   * `DeformationGradient` (SQCNI: the domain follows @f$ \boldsymbol{F}_n @f$), `None` (SNNI: the domain is only
   * translated with the particle), `RotationOnly` (the rotation @f$ \boldsymbol{R} @f$ of the polar decomposition)
   * and `RotationAndPrincipalStretch` (@f$ \boldsymbol{R} @f$ times the diagonal of
   * @f$ \boldsymbol{R}^\mathsf{T}\boldsymbol{F}_n @f$). Only for `DeformationGradient` does the smoothing domain
   * coincide with the physical domain of the particle.
   *
   * The volume, the smoothed gradients and the face loads are those of the section (ParticleSectionIntegrals): for
   * Section::Axisymmetric, the ring volume, the ring-weighted smoothed gradient with its area term, and two Gauss
   * points per face, both for the smoothing and for a load on the face. An axisymmetric particle needs the update type
   * `DeformationGradient` (smoothing domains that tile the body, for the integration constraint), refuses a vertex at
   * @f$ r < 0 @f$ (ParticleSectionIntegrals::checkParticle()), and supports the distributed load `PRESSURE` only. See
   * the page "Axisymmetric Particles" of the documentation.
   *
   * @tparam nDim      Spatial dimension (2: plane strain or axisymmetric, 3: 3D).
   * @tparam nVertices Number of vertices of the smoothing domain (4: Quad, 8: Hexa).
   * @tparam nNonlocal Number of nonlocal variables (1: a scalar nonlocal damage driver, 6: a nonlocal strain).
   * @tparam section   The section (Section::Cartesian: plane strain or 3D; Section::Axisymmetric: nDim = 2).
   */
  template < int nDim, int nVertices, int nNonlocal, Section section = Section::Cartesian >
  class GradientEnhancedFiniteStrainParticleSQCNI
    : public GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section > {

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
    using LagrangeCellType = LagrangeCell< nDim, nVertices >; ///< geometry of the smoothing domain
    /// the integrals over the smoothing domain
    using SectionIntegrals = ParticleSectionIntegrals< section, nDim, nVertices >;

    const SmoothingDomainUpdateType _smoothingVolumeUpdateType; ///< update type of the smoothing domain

    const Eigen::Matrix< double, nDim, nVertices > _vertexCoordinates_Undeformed; ///< undeformed vertex coordinates

    double*
      _vertexDisplacements_SmoothingDomain; ///< vertex displacements of the smoothing domain (in the state vector)

    double* _vertexDisplacements_Geometry;  ///< vertex displacements of the particle geometry (in the state vector)

    using ParentPointParticle = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >; ///< the point
                                                                                                  ///< particle base

    /// \brief The smoothing domain of the particle, built from its current vertex coordinates
    /// \return The cell of the smoothing domain in the intermediate configuration.
    LagrangeCellType smoothingDomain() const
    {
      Eigen::Matrix< double, nDim, nVertices > vertexCoordinates;
      getSmoothingDomainVertexCoordinates( vertexCoordinates.data() );
      return LagrangeCellType( vertexCoordinates.data(), nDim * nVertices );
    }

    /// \brief The geometry of the particle, built from its current vertex coordinates (getVertexCoordinates())
    /// \return The cell of the particle in the intermediate configuration; for the update type DeformationGradient it
    /// coincides with the smoothing domain.
    LagrangeCellType geometry() const
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

    /**
     * @brief Vertex coordinates of the smoothing domain: undeformed coordinates plus `smoothing vertex
     * displacements`.
     * @param[out] coordinates Coordinates, vertex by vertex (nDim x nVertices).
     */
    void getSmoothingDomainVertexCoordinates( double* coordinates ) const;

    /**
     * @brief The undeformed vertices mapped with a deformation tensor about the undeformed centroid and moved with the
     * displacement of the center: @f$ \boldsymbol{X}_c + \boldsymbol{A}\,(\boldsymbol{X}_v - \boldsymbol{X}_c) +
     * \boldsymbol{u}_c @f$.
     * @param[in] A The deformation tensor (@f$ \boldsymbol{F}_n @f$ for the geometry, see SmoothingDomainUpdateType
     * for the smoothing domain).
     * @return The vertex coordinates (nDim x nVertices).
     */
    Eigen::Matrix< double, nDim, nVertices > _vertexCoordinatesDeformedWith(
      const Eigen::Matrix< double, nDim, nDim >& A ) const;

  public:
    /**
     * @brief Vertex coordinates of the particle geometry: the undeformed vertices mapped with the deformation gradient
     * @f$ \boldsymbol{F}_n @f$ of the last accepted state about the undeformed centroid, plus the center displacement,
     * @f$ \boldsymbol{x}_v = \boldsymbol{X}_c + \boldsymbol{F}_n(\boldsymbol{X}_v - \boldsymbol{X}_c) +
     * \boldsymbol{u}_c @f$ (the undeformed vertices plus the state `vertex displacements`, updated at each accepted
     * increment). These are the points the particle moves with getInterpolationVectorLinearizedAboutCenter() (its own
     * points), the vertices supports and contacts act on. They coincide with the smoothing domain only for
     * `DeformationGradient`.
     * @param[out] coordinates Coordinates, vertex by vertex (nDim x nVertices).
     */
    void getVertexCoordinates( double* coordinates ) const override;

    /**
     * @brief Vertex coordinates for visualization: the smoothing domain (getSmoothingDomainVertexCoordinates()).
     * @param[out] coordinates Coordinates, vertex by vertex.
     */
    void getVisualizationVertexCoordinates( double* coordinates ) const override
    {
      getSmoothingDomainVertexCoordinates( coordinates );
    };

    /**
     * @brief Number of vertices.
     * @return nVertices.
     */
    int getNumberOfVertices() const override { return nVertices; };

    /**
     * @brief Shape of the particle.
     * @return Shape name of the undeformed smoothing domain cell.
     */
    std::string getParticleShape() const override { return _makeUndeformedCell().getCellShape(); }

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
     * @param[in] materialName              Name of a Marmot::Materials::MaterialGradientEnhancedFiniteStrain material.
     * @param[in] materialProperties        Material properties.
     * @param[in] sizeMaterialProperties    Number of material properties.
     * @param[in] approximation             Meshfree approximation used for the shape functions.
     * @param[in] smoothingVolumeUpdateType Update type of the smoothing domain.
     * @throws std::invalid_argument for an axisymmetric particle the section cannot integrate
     * (ParticleSectionIntegrals::checkParticle()).
     */
    GradientEnhancedFiniteStrainParticleSQCNI( int                             elementID,
                                               const double*                   nodeCoordinates,
                                               int                             nNodeCoordiantes,
                                               double                          volume,
                                               const std::string&              materialName,
                                               const double*                   materialProperties,
                                               int                             sizeMaterialProperties,
                                               const Approximation&            approximation,
                                               const SmoothingDomainUpdateType smoothingVolumeUpdateType );

    /**
     * @brief Evaluate the shape functions with the assigned kernel functions, with smoothed gradients.
     *
     * @f$ N_B @f$ is evaluated at the material point position, @f$ \partial N_B/\partial\boldsymbol{Y} @f$ as the
     * boundary integral over the current smoothing domain with one point per face (see the class description); the
     * test functions are set equal to the trial functions (before any VCI correction).
     */
    void updateShapeFunctions() override;

    /// \brief Get the smoothing volume of the particle (in its section)
    /// \return The smoothing volume of the particle
    /// \details The smoothing volume is computed as the volume of the element in an updated configuration
    ///         that is obtained by applying (parts of) the deformation gradient to the element in the undeformed
    ///         configuration. In general, this is not consistent with the actual deformation of the particle.
    ///
    virtual double getSmoothingVolume() const { return SectionIntegrals::volume( smoothingDomain() ); }

    /**
     * @brief Center coordinates: the material point position of the last accepted state.
     * @param[out] coordinates Coordinates (nDim values).
     */
    void getCenterCoordinates( double* coordinates ) const override { this->_mp.getCoordinatesAtCenter( coordinates ); }

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
     * @brief Volume of a cell spanned by given vertices, in the section (the ring volume of an axisymmetric particle).
     * @param[in] vertexCoordinates Vertex coordinates.
     * @return The volume.
     */
    virtual double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates ) const
    {
      return SectionIntegrals::volume( LagrangeCellType( vertexCoordinates.data(), nDim * nVertices ) );
    }

    /**
     * @brief Supported distributed loads, on the faces of the particle domain.
     * @return `PRESSURE` and `CWFCORRECTION`; `PRESSURE` only for an axisymmetric particle.
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      return SectionIntegrals::template supportedFaceLoads< ParentPointParticle::Pressure,
                                                            ParentPointParticle::CWFCorrection >();
    };

    /**
     * @brief Volume of the last accepted state, @f$ V_0\det\boldsymbol{F}_n @f$; for an axisymmetric particle the ring
     * volume of its geometry (geometry()).
     * @return The volume in the intermediate reference configuration (the current increment is not included).
     */
    virtual double getVolumeDeformed() const
    {
      if constexpr ( SectionIntegrals::isAxisymmetric )
        return SectionIntegrals::volume( geometry() );
      else
        return this->getVolumeUndeformed() * Fastor::determinant( this->_mp.dY_dX() );
    }

    /**
     * @brief Number of state variables: those of the point particle plus the vertex displacements of the geometry
     * and of the smoothing domain.
     * @return Required size of the state variable vector.
     */
    int getNumberOfRequiredStateVars() const override
    {
      return GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::getNumberOfRequiredStateVars() +
             2 * nDim * nVertices;
    };

    /**
     * @brief Assign the state variable vector; the last 2 x nDim x nVertices entries hold the vertex displacements of
     * the geometry, then those of the smoothing domain. Both are zero (the undeformed vertices) until the first
     * accepted increment, also before initializeYourself().
     * @param[in,out] stateVars  State variable vector.
     * @param[in]     nStateVars Its size.
     */
    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::assignStateVars( stateVars,
                                                                                         nStateVars -
                                                                                           2 * nDim * nVertices );
      _vertexDisplacements_Geometry        = stateVars + nStateVars - 2 * nDim * nVertices;
      _vertexDisplacements_SmoothingDomain = stateVars + nStateVars - nDim * nVertices;
    }

    /**
     * @brief Access a state variable; `vertex displacements` are those of the particle geometry,
     * `smoothing vertex displacements` those of the smoothing domain (as for the displacement particles).
     * @param[in] stateName Name of the state variable.
     * @param[in] qp        Evaluation point (unused).
     * @return View on the state variable.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qp ) const override
    {
      if ( stateName == "vertex displacements" )
        return { _vertexDisplacements_Geometry, nDim * nVertices };
      if ( stateName == "smoothing vertex displacements" )
        return { _vertexDisplacements_SmoothingDomain, nDim * nVertices };
      return GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::getStateView( stateName, qp );
    }

    /**
     * @brief Distributed loads on a face of the particle.
     *
     * Both loads use the face quadrature in the intermediate reference configuration
     * (faceQuadratureIntermediate()): on each of its points @f$ \boldsymbol{Y}_f @f$, with its area vector
     * @f$ \boldsymbol{N}\,dA_Y @f$, the test functions @f$ T_A(\boldsymbol{Y}_f) @f$ evaluated by the meshfree
     * approximation. The loads are subtracted from @p fExt (sign convention of the internal force), with the consistent
     * tangent with respect to the displacement increment:
     *
     * - `PRESSURE` (follower pressure @f$ p @f$ = `load[0]`):
     *   @f$ \boldsymbol{f} = -\Delta J\,\Delta\boldsymbol{F}^{-\mathsf{T}}\,p\,\boldsymbol{N}\,dA_Y @f$,
     *   @f$ f_{Ai} \mathrel{-}= T_A(\boldsymbol{Y}_f)\,f_i @f$; as for the elements, a positive pressure pushes on
     *   the face (against the outward normal); see ParticleSectionIntegrals::addPressure(). For an axisymmetric
     *   particle, the face carries two Gauss points, the rule of the smoothing, and the pressure force of each is that
     *   of the half plane times the ring of the point after the increment (Hoop::ringOfFacePoint()), with its
     *   derivative in the tangent.
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
     * @throws std::invalid_argument for an unknown load type, and for `CWFCORRECTION` of an axisymmetric particle.
     */
    void computeDistributedLoad( int           type,
                                 int           surfaceID,
                                 const double* load,
                                 double*       fExt,
                                 double*       dExt_dQ,
                                 double        timeNew,
                                 double        dT ) const override;

    /**
     * @brief The quadrature of a face in the intermediate reference configuration, for the loads on it.
     *
     * For `DeformationGradient` (the only update of an axisymmetric particle) it is the rule of the section on the
     * face of the current smoothing domain, which coincides with the particle
     * (ParticleSectionIntegrals::faceQuadrature()): a face load is integrated with the rule of the smoothing. For the
     * other update types (plane strain and 3D) it is one point: the undeformed surface vector mapped by Nanson's
     * formula,
     * @f$ \boldsymbol{N}\,dA_Y = J_n\,\boldsymbol{F}_n^{-\mathsf{T}}\,\boldsymbol{N}\,dA_0 @f$, at the undeformed
     * face center translated by the particle displacement.
     *
     * @param[in] boundaryFaceID Face id (1-based).
     * @return The points and their share in the area vector of the face.
     */
    typename SectionIntegrals::FaceQuadrature faceQuadratureIntermediate( int boundaryFaceID ) const;

    /**
     * @brief Add the boundary term @f$ T_A(\boldsymbol{Y}_f)\,P_C(\boldsymbol{Y}_f)\,N_i\,dA_Y @f$ of the VCI
     * integration constraint for a face of the particle.
     *
     * Overrides the point version: the points and their surface vectors come from faceQuadratureIntermediate(), and
     * @f$ T_A @f$ and @f$ P_C @f$ are evaluated at the points.
     *
     * @param[in,out] R_AiC_RowMajor        Constraint residual @f$ R_{AiC} @f$, row-major.
     * @param[in]     boundarySurfaceVector Unused.
     * @param[in]     boundaryFaceID        Face id (1-based).
     */
    void vci_compute_Test_P_BoundaryIntegral( double*       R_AiC_RowMajor,
                                              const double* boundarySurfaceVector,
                                              int           boundaryFaceID ) override
    {
      const auto face = faceQuadratureIntermediate( boundaryFaceID );
      for ( int g = 0; g < SectionIntegrals::nPointsPerFace; g++ ) {
        const CoordinatesSized Y_N       = face.points.col( g );
        const CoordinatesSized N_dAY     = SectionIntegrals::weight( Y_N ) * face.n_dA;
        Eigen::MatrixXd        TBoundary = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

        ParentPointParticle::_meshfreeApproximation
          .computeShapeFunctions( Y_N.data(), ParentPointParticle::_assignedKernelFunctions, TBoundary.data() );

        // get P for the exact integration location at the boundary
        auto PBoundary = ParentPointParticle::_P;
        Math::computeMonomialBasis( ParentPointParticle::_vciOrder, Y_N, PBoundary );

        for ( int A = 0; A < ParentPointParticle::_nNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            for ( int C = 0; C < ParentPointParticle::_nVCIConstraints; C++ )
              R_AiC_RowMajor[A * ( nDim * ParentPointParticle::_nVCIConstraints ) +
                             i * ParentPointParticle::_nVCIConstraints + C] += TBoundary( A ) * PBoundary( C ) *
                                                                               N_dAY[i];
      }
    };

    /**
     * @brief Coordinates of the evaluation points: the quadrature points of the faces of the smoothing domain (the face
     * centers, or two Gauss points per face of an axisymmetric particle,
     * ParticleSectionIntegrals::faceQuadraturePoints()).
     * @param[out] coordinates Coordinates, face by face (nDim x number of evaluation points).
     */
    void getEvaluationCoordinates( double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > points( coordinates,
                                                                          nDim,
                                                                          getNumberOfEvaluationPoints() );
      points = SectionIntegrals::faceQuadraturePoints( smoothingDomain() );
    }

    /**
     * @brief Center of a face of the particle geometry (geometry()).
     * @param[in]  faceID      Face id (1-based).
     * @param[out] coordinates Face center (nDim values).
     */
    void getFaceCoordinates( int faceID, double* coordinates ) const override
    {
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > > faceCenter( coordinates );
      faceCenter = geometry().getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Number of evaluation points.
     * @return The number of face quadrature points of the smoothing domain (whose faces are those of the undeformed
     * cell).
     */
    int getNumberOfEvaluationPoints() const override
    {
      return SectionIntegrals::numberOfFaceQuadraturePoints( _makeUndeformedCell() );
    };

  private:
    /// @brief Update the vertex displacements of the smoothing domain from @f$ \boldsymbol{F}_n @f$ and the particle
    /// displacement, according to the SmoothingDomainUpdateType, and those of the geometry with @f$ \boldsymbol{F}_n
    /// @f$.
    void _updateVertexDisplacementsFromMaterialPointDeformation();

    /// @brief Update the smoothing domain and the geometry, then set the center to the accepted material point
    /// position.
    void updateParticlePositionToReferenceIntermediate() override
    {
      _updateVertexDisplacementsFromMaterialPointDeformation();

      getCenterCoordinates( ParentPointParticle::_centerReferenceIntermediate.data() );
    };
  };

  template < int nDim, int nVertices, int nNonlocal, Section section >
  GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::
    GradientEnhancedFiniteStrainParticleSQCNI( int                elementID,
                                               const double*      vertexCoordinates,
                                               int                nVertexCoordinates,
                                               double             volume,
                                               const std::string& materialName,
                                               const double*      materialProperties,
                                               int                sizeMaterialProperties,
                                               const Marmot::Discretizations::Meshfree::Approximation& approximation,
                                               const SmoothingDomainUpdateType smoothingVolumeUpdateType )
    : GradientEnhancedFiniteStrainParticle< nDim,
                                            nNonlocal,
                                            section >( elementID,
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
    SectionIntegrals::checkParticle( _makeUndeformedCell(),
                                     smoothingVolumeUpdateType == DeformationGradient,
                                     elementID );
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::computeDistributedLoad(
    int           type,
    int           boundaryFaceID,
    const double* load_,
    double*       fExt,
    double*       dFExt_ddQ,
    double        timeNew,
    double        dT ) const
  {

    if ( SectionIntegrals::isAxisymmetric &&
         type == GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::CWFCorrection )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": CWFCORRECTION is not available for an axisymmetric particle" );

    switch ( type ) {

    case GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::Pressure: {
      SectionIntegrals::template addPressure< nDim + nNonlocal >( faceQuadratureIntermediate( boundaryFaceID ),
                                                                  load_[0],
                                                                  ParentPointParticle::_mp,
                                                                  ParentPointParticle::_N,
                                                                  ParentPointParticle::_dN_dY,
                                                                  ParentPointParticle::_meshfreeApproximation,
                                                                  ParentPointParticle::_assignedKernelFunctions,
                                                                  fExt,
                                                                  dFExt_ddQ );
      break;
    }
    case GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::CWFCorrection: {
      const auto&   _mp           = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp;
      const auto&   _nNodes       = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_nNodes;
      constexpr int nodeBlockSize = nDim + nNonlocal;

      Eigen::Map< Eigen::VectorXd > P( fExt, _nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, _nNodes * nodeBlockSize, _nNodes * nodeBlockSize );

      using namespace Fastor;

      const auto face = faceQuadratureIntermediate( boundaryFaceID );
      for ( int g = 0; g < SectionIntegrals::nPointsPerFace; g++ ) {
        const CoordinatesSized Y_N = face.points.col( g );
        const TensorD          N_dAY( CoordinatesSized( SectionIntegrals::weight( Y_N ) * face.n_dA ).data() );

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

        std::array< TensorD, nNonlocal > dTraction_dN;
        for ( int kNl = 0; kNl < nNonlocal; kNl++ )
          dTraction_dN[kNl] = t.dS_dN[kNl] % v;

        // load[0] selects the corrected components as a bit mask (1 = x, 2 = y, 4 = z; 0 or no load = all): on a face
        // where only the normal displacement is prescribed (symmetry plane, frictionless platen) only that component of
        // the traction is a reaction, the tangential one is a natural (zero) condition
        const int mask = load_ ? static_cast< int >( std::lround( load_[0] ) ) : 0;
        TensorD   sel;
        for ( int i = 0; i < nDim; ++i )
          sel( i ) = ( mask == 0 || ( mask >> i ) & 1 ) ? 1.0 : 0.0;
        for ( int i = 0; i < nDim; ++i ) {
          traction( i ) *= sel( i );
          for ( int kNl = 0; kNl < nNonlocal; kNl++ )
            dTraction_dN[kNl]( i ) *= sel( i );
          for ( int m = 0; m < nDim; ++m )
            for ( int M = 0; M < nDim; ++M )
              dTraction_dDeltaF( i, m, M ) *= sel( i );
        }

        Eigen::RowVectorXd testBoundary = Eigen::RowVectorXd::Zero( ParentPointParticle::_nNodes );

        ParentPointParticle::_meshfreeApproximation
          .computeShapeFunctions( Y_N.data(), ParentPointParticle::_assignedKernelFunctions, testBoundary.data() );

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
              GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_dN_dY.col( B ).data() );
            const double N_B = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_N( B );

            Tensor< double, nDim, nDim > k_UU( 0.0 );
            for ( int i = 0; i < nDim; ++i )
              for ( int m = 0; m < nDim; ++m )
                for ( int M = 0; M < nDim; ++M )
                  k_UU( i, m ) += testBoundary( A ) * dTraction_dDeltaF( i, m, M ) * dN_B_dY( M );

            {
              using namespace Eigen;
              K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< const Matrix< double, nDim, nDim > >(
                torowmajor( k_UU ).data() );
              for ( int kNl = 0; kNl < nNonlocal; kNl++ ) {
                const TensorD k_UN = testBoundary( A ) * N_B * dTraction_dN[kNl];
                K.template block< nDim, 1 >( idxA_u,
                                             idxB_n + kNl ) -= Map< const Matrix< double, nDim, 1 > >( k_UN.data() );
              }
            }
          }
        }
      }
      break;
    }
    default: {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  typename GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::SectionIntegrals::
    FaceQuadrature
    GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::faceQuadratureIntermediate(
      int boundaryFaceID ) const
  {
    if constexpr ( !SectionIntegrals::isAxisymmetric )
      if ( _smoothingVolumeUpdateType != DeformationGradient ) {
        // The deformed smoothing domain is not consistent with the physical domain of the particle, so the face is
        // not taken from it: the undeformed boundary vector is mapped by Nanson's formula with F_n, and placed at the
        // undeformed face center translated by the particle displacement.
        const auto                             cellUndeformed = _makeUndeformedCell();
        const Eigen::Matrix< double, nDim, 1 > N_dA0_eigen    = cellUndeformed.boundarySurfaceVector( boundaryFaceID );
        const Eigen::Matrix< double, nDim, 1 > Y0 = cellUndeformed.getFaceCenterCoordinates( boundaryFaceID );

        const auto FY = this->_mp.dY_dX();

        using namespace Fastor;
        const TensorD N_dAY = determinant( FY ) * transpose( inverse( FY ) ) % TensorD( N_dA0_eigen.data() );

        TensorD Y;
        this->_mp.getCenterDisplacement( Y.data() );
        Y += TensorD( Y0.data() );

        return { CoordinatesSized( Y.data() ), CoordinatesSized( N_dAY.data() ) };
      }

    // the full SQCNI: the deformed smoothing domain coincides with the physical domain of the particle
    return SectionIntegrals::faceQuadrature( smoothingDomain(), boundaryFaceID );
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::
    getSmoothingDomainVertexCoordinates( double* coordinates ) const
  {

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > coordinatesDeformed( coordinates );
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacementsSmoothingDomain(
      _vertexDisplacements_SmoothingDomain );

    coordinatesDeformed = _vertexCoordinates_Undeformed + vertexDisplacementsSmoothingDomain;
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::
    _updateVertexDisplacementsFromMaterialPointDeformation()
  {

    Eigen::Matrix< double, nDim, nDim > F = Eigen::Matrix< double, nDim, nDim >::Identity();

    switch ( _smoothingVolumeUpdateType ) {

    case SmoothingDomainUpdateType::None: F.setIdentity(); break;

    case SmoothingDomainUpdateType::DeformationGradient: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp
                                                  .dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );
      break;
    }

    case SmoothingDomainUpdateType::RotationOnly: {

      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp
                                                  .dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      const Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::ComputeFullU | Eigen::ComputeFullV > svd( F );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      F = R;
      break;
    }

    case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
      Fastor::Tensor< double, nDim, nDim > F_ = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp
                                                  .dY_dX();
      F = Eigen::Map< Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >( F_.data() );

      const Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::ComputeFullU | Eigen::ComputeFullV > svd( F );

      Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();

      Eigen::Matrix< double, nDim, nDim > U_ = Eigen::Matrix< double, nDim, nDim >::Identity();
      U_.diagonal()                          = ( R.transpose() * F ).diagonal();

      F = R * U_;

      break;
    }
    }

    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacements( _vertexDisplacements_SmoothingDomain );
    vertexDisplacements = _vertexCoordinatesDeformedWith( F ) - _vertexCoordinates_Undeformed;

    const Fastor::Tensor< double, nDim, nDim >
      F_n = GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp.dY_dX();
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > > vertexDisplacementsGeometry( _vertexDisplacements_Geometry );
    vertexDisplacementsGeometry = _vertexCoordinatesDeformedWith(
                                    Eigen::Map< const Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor > >(
                                      F_n.data() ) ) -
                                  _vertexCoordinates_Undeformed;
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  Eigen::Matrix< double, nDim, nVertices > GradientEnhancedFiniteStrainParticleSQCNI<
    nDim,
    nVertices,
    nNonlocal,
    section >::_vertexCoordinatesDeformedWith( const Eigen::Matrix< double, nDim, nDim >& A ) const
  {
    Eigen::Matrix< double, nDim, nVertices > coordinatesRelToCenter0 = _vertexCoordinates_Undeformed;
    for ( int i = 0; i < nDim; i++ ) {
      coordinatesRelToCenter0.row( i ).array() -= this->_centerCoordinatesUndeformed( i );
    }

    Eigen::Matrix< double, nDim, nVertices > vertexCoordinatesDeformed = A * coordinatesRelToCenter0;

    Eigen::Matrix< double, nDim, 1 > _mp_displacement;
    GradientEnhancedFiniteStrainParticle< nDim, nNonlocal, section >::_mp.getCenterDisplacement(
      _mp_displacement.data() );

    for ( int i = 0; i < nDim; i++ ) {
      vertexCoordinatesDeformed.row( i ).array() += this->_centerCoordinatesUndeformed( i ) + _mp_displacement( i );
    }

    return vertexCoordinatesDeformed;
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::getVertexCoordinates(
    double* coordinates ) const
  {
    Eigen::Map< Eigen::Matrix< double, nDim, nVertices > >       coordinatesDeformed( coordinates );
    Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > > vertexDisplacements( _vertexDisplacements_Geometry );
    coordinatesDeformed = _vertexCoordinates_Undeformed + vertexDisplacements;
  }

  template < int nDim, int nVertices, int nNonlocal, Section section >
  void GradientEnhancedFiniteStrainParticleSQCNI< nDim, nVertices, nNonlocal, section >::updateShapeFunctions()
  {
    Eigen::Matrix< double, nDim, 1 > coords;
    ParentPointParticle::_mp.getCoordinatesAtCenter( coords.data() );

    ParentPointParticle::_N = Eigen::MatrixXd::Zero( 1, ParentPointParticle::_nNodes );

    ParentPointParticle::_meshfreeApproximation.computeShapeFunctions( coords.data(),
                                                                       ParentPointParticle::_assignedKernelFunctions,
                                                                       ParentPointParticle::_N.data() );

    SectionIntegrals::smoothShapeFunctionGradients( smoothingDomain(),
                                                    ParentPointParticle::_meshfreeApproximation,
                                                    ParentPointParticle::_assignedKernelFunctions,
                                                    ParentPointParticle::_N,
                                                    ParentPointParticle::_dN_dY,
                                                    nullptr ); // no second derivatives (not NSNI)

    ParentPointParticle::_T     = ParentPointParticle::_N;
    ParentPointParticle::_dT_dY = ParentPointParticle::_dN_dY;
  }

} // namespace Marmot::Discretizations::Meshfree
