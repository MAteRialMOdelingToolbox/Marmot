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

#include "Marmot/MarmotMeshfreeQuadHexCell.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <vector> // Explicitly include vector for clarity

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::ParticleDomain
   * @brief The geometry of a cell-shaped particle, independent of physics or meshfree approximation.
   * @details A particle domain is a Lagrangian cell (MarmotLagrangeCell, e.g., a 4-node quadrilateral or an 8-node
   * hexahedron) that is kept in three versions:
   *  - the **undeformed geometry** with the vertices @f$ \boldsymbol{X}_v @f$ and the centroid
   *    @f$ \boldsymbol{X}_c @f$ given at construction,
   *  - the **deformed geometry** (the geometry in the reference intermediate configuration, i.e., the last accepted
   *    configuration), used for the particle position, the face centers, the boundary surface vectors of
   *    distributed loads and the second moments of the NSNI stabilization,
   *  - the **smoothing domain**, over whose boundary the smoothed shape function gradients are integrated (see
   *    GenericSDIParticle and the SQCNI particles).
   *
   * Both deformed versions are obtained from the undeformed geometry by a homogeneous deformation about the
   * centroid plus the displacement of the center (see acceptStateAndPosition()),
   * @f[
   *   \boldsymbol{x}_v = \boldsymbol{X}_c + \boldsymbol{F}\,(\boldsymbol{X}_v - \boldsymbol{X}_c)
   *   + \boldsymbol{u}_c ,
   * @f]
   * where @f$ \boldsymbol{X}_c @f$ is the centroid of the particle; for the subdomains of uniformSubdivided() it is
   * the centroid of the parent domain, so that the subdomains remain a tiling of the deformed parent.
   * with the total deformation gradient @f$ \boldsymbol{F} @f$ of the particle for the deformed geometry, and a
   * tensor @f$ \boldsymbol{F}_s @f$ derived from it according to SmoothingDomainUpdateType for the smoothing domain.
   * Face IDs are 1-based and follow the Abaqus convention of MarmotLagrangeCell.
   *
   * @tparam nDim The number of dimensions (2 or 3).
   * @tparam nVertices The number of vertices of the cell (4 for a quadrilateral, 8 for a hexahedron).
   */
  template < int nDim, int nVertices >
  class ParticleDomain {

  public:
    // Type aliases for improved readability
    using CoordinatesSized       = Eigen::Matrix< double, nDim, 1 >; ///< Alias for an Eigen vector storing coordinates.
    using VertexCoordinatesSized = Eigen::Matrix< double, nDim, nVertices >; ///< Alias for Eigen matrix storing vertex
                                                                             ///< coordinates.
    using DeformationGradientSized = Eigen::Matrix< double, nDim, nDim >;    ///< Alias for Eigen matrix storing
                                                                             ///< deformation gradient.
    using LagrangeCellType = MarmotLagrangeCell< nDim, nVertices >; ///< Alias for the underlying Lagrange cell type.

    /**
     * @brief Defines how the smoothing domain follows the deformation of the particle.
     * @details The smoothing domain is mapped with the tensor @f$ \boldsymbol{F}_s @f$ computed by
     * _computeSmoothingDomainDeformationTensorTotal(); in all cases it is translated by the center displacement.
     * The registered particle names encode the choice (e.g., SQCNI: DeformationGradient, SNNI: None, suffix _R or
     * prefix R-: RotationOnly, suffix _RU or prefix RS-: RotationAndPrincipalStretch).
     */
    enum SmoothingDomainUpdateType {
      None,                       ///< @f$ \boldsymbol{F}_s = \boldsymbol{I} @f$: the smoothing domain keeps its
                                  ///< undeformed shape and is only translated.
      DeformationGradient,        ///< @f$ \boldsymbol{F}_s = \boldsymbol{F} @f$: the smoothing domain deforms with
                                  ///< the particle (it coincides with the deformed geometry).
      RotationOnly,               ///< @f$ \boldsymbol{F}_s = \boldsymbol{R} @f$, the rotation of the polar
                                  ///< decomposition @f$ \boldsymbol{F} = \boldsymbol{R}\boldsymbol{U} @f$ (from an
                                  ///< SVD).
      RotationAndPrincipalStretch ///< @f$ \boldsymbol{F}_s = \boldsymbol{R}\,\mathrm{diag}(\boldsymbol{R}^T
                                  ///< \boldsymbol{F}) @f$: the rotation and the diagonal entries of
                                  ///< @f$ \boldsymbol{U} @f$ in the global basis.
    };

    /**
     * @brief Constructs a new ParticleDomain object; all three geometries start as the undeformed geometry.
     * @param[in] vertexCoordinates Pointer to an array of vertex coordinates (nDim * nVertices, vertex by vertex)
     * in the undeformed configuration.
     * @param[in] nVertexCoordinates The total number of coordinate values (nDim * nVertices).
     * @param[in] smoothingVolumeUpdateType The strategy for updating the smoothing domain.
     */
    ParticleDomain( const double*                   vertexCoordinates,
                    int                             nVertexCoordinates,
                    const SmoothingDomainUpdateType smoothingVolumeUpdateType );

    const SmoothingDomainUpdateType smoothingVolumeUpdateType; ///< Type of update of the smoothing domain.

    /**
     * @brief Gets the coordinates of the particle's vertices in the intermediate (deformed) configuration.
     * @return A const reference to an Eigen matrix containing the vertex coordinates.
     */
    const VertexCoordinatesSized& getGeometryDeformedVertexCoordinates() const
    {
      return _cellForGeometryDeformed.nodes();
    }

    /**
     * @brief Gets the coordinates of the particle's smoothing domain vertices.
     * @return A const reference to an Eigen matrix containing the vertex coordinates.
     */
    const VertexCoordinatesSized& getSmoothingVertexCoordinates() const { return _cellForSmoothing.nodes(); }

    /**
     * @brief Gets the number of vertices defining the particle.
     * @return The number of vertices.
     */
    int getNumberOfVertices() const { return nVertices; }

    /**
     * @brief Gets the shape of the particle (e.g., "quad", "hex").
     * @return A string representing the particle's shape.
     */
    std::string getParticleShape() const { return _cellForGeometryUndeformed.getCellShape(); }

    /**
     * @brief Gets the coordinates of a specific face's center in the intermediate (deformed) configuration.
     *        This refers to the faces of the *deformed geometry*.
     * @param[in] faceID The ID of the face (1-based index).
     * @return An Eigen vector representing the face center coordinates.
     */
    CoordinatesSized getFaceCenterCoordinates( int faceID ) const
    {
      return _cellForGeometryDeformed.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Gets the coordinates of a specific evaluation point (face center of the smoothing domain).
     * @param[in] faceID The ID of the face (1-based index).
     * @return An Eigen vector representing the evaluation point coordinate.
     */
    CoordinatesSized getSmoothingDomainFaceCenterCoordinates( int faceID ) const
    {
      return _cellForSmoothing.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Gets the number of faces of the cell (the same for all three geometries); this is also the number of
     *        evaluation points of the smoothing boundary integral.
     * @return The number of faces.
     */
    int getNumberOfFaces() const { return _cellForGeometryDeformed.getNumberOfFaces(); }

    /**
     * @brief Get the volume of the smoothing domain in its current state.
     * @return The smoothing volume of the particle.
     */
    double getSmoothingVolume() const { return _cellForSmoothing.volume(); }

    /**
     * @brief Gets the center coordinates of the particle in the intermediate (deformed) configuration.
     * @return An Eigen vector representing the center coordinates.
     */
    CoordinatesSized getCenterCoordinates() const { return _cellForGeometryDeformed.centroid(); }

    /**
     * @brief Gets the undeformed volume of the particle.
     * @return The undeformed volume.
     */
    double getVolumeUndeformed() const { return _cellForGeometryUndeformed.volume(); }

    /**
     * @brief Provides a const reference to the vertex displacements of the smoothing domain (smoothing domain
     *        vertices minus undeformed vertices).
     * @return A const reference to an Eigen matrix containing the vertex displacements.
     */
    const VertexCoordinatesSized& getSmoothingDomainVertexDisplacements() const
    {
      return _vertex_displacements_smoothingDomain;
    }

    /**
     * @brief Provides a const reference to the vertex displacements of the geometry (deformed minus undeformed
     *        vertices).
     * @return A const reference to an Eigen matrix containing the vertex displacements.
     */
    const VertexCoordinatesSized& getGeometryDeformedVertexDisplacements() const
    {
      return _vertex_displacements_geometry;
    }

    /**
     * @brief Returns the boundary surface vector for a given face ID, from the *deformed geometry*.
     *        This is typically used for distributed loads.
     * @details The vector is the outward normal scaled by the face area (edge length in 2D), @f$ \boldsymbol{n}\,dA
     * @f$.
     * @param[in] faceID The ID of the face (1-based index).
     * @return An Eigen vector representing the boundary surface vector.
     */
    CoordinatesSized getFaceBoundaryVector( int faceID ) const
    {
      return _cellForGeometryDeformed.boundarySurfaceVector( faceID );
    }

    /**
     * @brief Returns the indices of sub-cells (of uniformSubdivided()) that lie on a given parent face.
     * @param[in] parentFaceId The ID of the parent face (1-based).
     * @return A vector of integers representing the sub-cell indices.
     */
    inline std::vector< int > getSubCellIndicesOnParentFace( int parentFaceId ) const
    {
      return _cellForGeometryUndeformed.getSubCellIndicesOnParentFace( parentFaceId );
    }

    /**
     * @brief Returns the boundary surface vector for a given face ID, from the *smoothing domain*.
     *        This is typically used for computing smoothed shape function gradients.
     * @param[in] faceID The ID of the face (1-based index).
     * @return An Eigen vector representing the boundary surface vector.
     */
    CoordinatesSized getSmoothingBoundarySurfaceVector( int faceID ) const
    {
      return _cellForSmoothing.boundarySurfaceVector( faceID );
    }

    /**
     * @brief Returns the second moments of area/volume for the deformed geometry,
     *        @f$ \int (\boldsymbol{x} - \boldsymbol{x}_c) \otimes (\boldsymbol{x} - \boldsymbol{x}_c)\, dV @f$ about
     * the centroid (Gauss quadrature on the cell).
     * @return An Eigen matrix representing the second moments.
     */
    DeformationGradientSized getGeometrySecondMoments() const { return _cellForGeometryDeformed.secondMoments(); }

    /**
     * @brief Updates the particle's position and volume to the reference intermediate configuration.
     * @details The deformed geometry and the smoothing domain are rebuilt from the undeformed geometry,
     *          @f$ \boldsymbol{x}_v = \boldsymbol{X}_c + \boldsymbol{F}(\boldsymbol{X}_v - \boldsymbol{X}_c) +
     *          \boldsymbol{u}_c @f$ (@f$ \boldsymbol{X}_c @f$: the centroid, of the parent domain for a subdomain),
     *          with @f$ \boldsymbol{F} @f$ = @p F_physics for the geometry and
     *          @f$ \boldsymbol{F}_s @f$ (see SmoothingDomainUpdateType) for the smoothing domain. The vertex
     *          displacements are updated accordingly.
     * @param[in] F_physics The total deformation gradient of the particle (with respect to the undeformed
     * configuration).
     * @param[in] centerDisplacement The total displacement of the particle's center.
     */
    void acceptStateAndPosition( const DeformationGradientSized& F_physics, const CoordinatesSized& centerDisplacement )
    {
      _cellForGeometryDeformed.updateVertexCoordinates( _cellForGeometryUndeformed.nodes() );
      _cellForGeometryDeformed.applyDeformationGradient( F_physics, _deformationCenter );
      _cellForGeometryDeformed.applyUniformDisplacement( centerDisplacement );

      const auto FSmoothing = _computeSmoothingDomainDeformationTensorTotal( F_physics );

      _vertex_displacements_geometry = _cellForGeometryDeformed.nodes() - _cellForGeometryUndeformed.nodes();

      _cellForSmoothing.updateVertexCoordinates( _cellForGeometryUndeformed.nodes() );
      _cellForSmoothing.applyDeformationGradient( FSmoothing, _deformationCenter );
      _cellForSmoothing.applyUniformDisplacement( centerDisplacement );

      _vertex_displacements_smoothingDomain = _cellForSmoothing.nodes() - _cellForGeometryUndeformed.nodes();
    }

    /**
     * @brief Uniformly subdivides the particle domain into smaller domains.
     * @details This method creates a vector of new `ParticleDomain` instances, each representing
     *          a uniformly subdivided portion of the original undeformed particle domain (one level of bisection
     *          in each direction, i.e., 4 quadrilaterals or 8 hexahedra), with the same smoothing domain update
     *          type. The subdomains are deformed about the centroid of this domain (not their own), so that the
     *          deformed subdomains tile the deformed domain.
     * @return A vector of `ParticleDomain` instances representing the subdivided particles.
     */
    std::vector< ParticleDomain > uniformSubdivided() const
    {
      std::vector< ParticleDomain > subdividedDomains;

      auto subdividedCells = _cellForGeometryUndeformed.uniformSubdivided();

      for ( const auto& cell : subdividedCells ) {
        ParticleDomain newDomain( cell.nodes().data(), nDim * nVertices, smoothingVolumeUpdateType );
        newDomain._deformationCenter = _deformationCenter;
        subdividedDomains.push_back( newDomain );
      }

      return subdividedDomains;
    }

    /**
     * @brief Computes the centroid coordinates from a given set of vertex coordinates.
     * @param[in] vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return An Eigen vector representing the centroid coordinates.
     */
    static CoordinatesSized getCenterFromVertices( const VertexCoordinatesSized& vertexCoordinates )
    {
      LagrangeCellType _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.centroid();
    }

    /**
     * @brief Computes the volume of a cell defined by a given set of vertex coordinates.
     * @param[in] vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return The computed volume.
     */
    static double getVolumeFromVertices( const VertexCoordinatesSized& vertexCoordinates )
    {
      LagrangeCellType _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.volume();
    }

  protected:
    LagrangeCellType _cellForGeometryUndeformed; ///< Cell representing the undeformed geometry.
    LagrangeCellType _cellForGeometryDeformed;   ///< Cell representing the intermediate (deformed) geometry.
    LagrangeCellType _cellForSmoothing;          ///< Cell representing the smoothing domain.

    VertexCoordinatesSized
      _vertex_displacements_smoothingDomain;               ///< Displacements of vertices in the smoothing domain.
    VertexCoordinatesSized _vertex_displacements_geometry; ///< Displacements of vertices of the geometry.
    CoordinatesSized       _deformationCenter; ///< Fixed point @f$ \boldsymbol{X}_c @f$ of the homogeneous deformation.

    /**
     * @brief Computes the total deformation tensor @f$ \boldsymbol{F}_s @f$ for the smoothing domain based on the
     *        configured update type (see SmoothingDomainUpdateType).
     * @param[in] F_physics The total deformation gradient of the particle.
     * @return An Eigen matrix representing the deformation tensor for the smoothing domain.
     */
    DeformationGradientSized _computeSmoothingDomainDeformationTensorTotal(
      const DeformationGradientSized& F_physics ) const
    {
      DeformationGradientSized F_smoothing;

      switch ( smoothingVolumeUpdateType ) {
      case SmoothingDomainUpdateType::None: F_smoothing.setIdentity(); break;
      case SmoothingDomainUpdateType::DeformationGradient: F_smoothing = F_physics; break;
      case SmoothingDomainUpdateType::RotationOnly: {
        Eigen::JacobiSVD< Eigen::MatrixXd > svd;
        svd.compute( F_physics, Eigen::ComputeFullU | Eigen::ComputeFullV );
        F_smoothing = svd.matrixU() * svd.matrixV().transpose();
        break;
      }
      case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
        Eigen::JacobiSVD< Eigen::MatrixXd > svd;
        svd.compute( F_physics, Eigen::ComputeFullU | Eigen::ComputeFullV );
        DeformationGradientSized R  = svd.matrixU() * svd.matrixV().transpose();
        DeformationGradientSized U_ = DeformationGradientSized::Identity();
        U_.diagonal()               = ( R.transpose() * F_physics ).diagonal();
        F_smoothing                 = R * U_;
        break;
      }
      }
      return F_smoothing;
    }
  };

  /*
   * (Definition; documented at the declaration.)
   * @brief Constructor for ParticleDomain.
   * @tparam nDim The number of dimensions.
   * @tparam nVertices The number of vertices.
   * @param[in] vertexCoordinates Pointer to an array of vertex coordinates (nDim * nVertices) in the undeformed
   * configuration.
   * @param[in] nVertexCoordinates The total number of coordinate values (nDim * nVertices).
   * @param[in] smoothingVolumeUpdateType The strategy for updating the smoothing domain.
   */
  template < int nDim, int nVertices >
  ParticleDomain< nDim, nVertices >::ParticleDomain( const double*                   vertexCoordinates,
                                                     int                             nVertexCoordinates,
                                                     const SmoothingDomainUpdateType smoothingVolumeUpdateType )

    : smoothingVolumeUpdateType( smoothingVolumeUpdateType ),
      _cellForGeometryUndeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForGeometryDeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForSmoothing( vertexCoordinates, nVertexCoordinates ),
      _vertex_displacements_smoothingDomain( VertexCoordinatesSized::Zero() ),
      _vertex_displacements_geometry( VertexCoordinatesSized::Zero() ),
      _deformationCenter( _cellForGeometryUndeformed.centroid() )
  {
  }

} // namespace Marmot::Meshfree
