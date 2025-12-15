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

namespace Marmot::Meshfree {

  /**
   * @brief A generic class for cell-based particle geometry, independent of physics or meshfree approximation.
   *
   * This class provides the geometric representation and operations for particles
   * defined by a cell (e.g., quadrilateral, hexahedral). It manages undeformed,
   * intermediate, and smoothing domain geometries.
   *
   * @tparam nDim The number of dimensions (e.g., 2 for 2D, 3 for 3D).
   * @tparam nVertices The number of vertices defining the particle's geometry.
   */
  template < int nDim, int nVertices >
  class ParticleDomain {

  public:
    // Enums and Typedefs
    using CoordinatesSized = Eigen::Matrix< double, nDim, 1 >;

    /**
     * @brief Defines how the smoothing domain's volume is updated.
     */
    enum SmoothingDomainUpdateType {
      None,                       ///< No update to the smoothing domain's deformation tensor (identity).
      DeformationGradient,        ///< Update using the total deformation gradient F.
      RotationOnly,               ///< Update using only the rotation part R from F = RU.
      RotationAndPrincipalStretch ///< Update using the rotation R and principal stretches from U.
    };

    // Constructor
    /**
     * @brief Constructs a new ParticleDomain object.
     *
     * @param vertexCoordinates Pointer to an array of vertex coordinates (nDim * nVertices).
     * @param nVertexCoordinates The total number of coordinate values (nDim * nVertices).
     * @param smoothingVolumeUpdateType The strategy for updating the smoothing domain's volume.
     */
    ParticleDomain( const double*                      vertexCoordinates,
                    int                                nVertexCoordinates,
                    const SmoothingDomainUpdateType    smoothingVolumeUpdateType );

    // Public Member Variables
    SmoothingDomainUpdateType smoothingVolumeUpdateType; ///< Type of update for the smoothing volume.

    // Public Methods (Getters)
    /**
     * @brief Gets the coordinates of the particle's vertices in the intermediate configuration.
     * @return A const reference to an Eigen matrix containing the vertex coordinates.
     */
    const Eigen::Matrix<double, nDim, nVertices>& getGeometryDeformedVertexCoordinates() const
    {
      return _cellForGeometryDeformed.nodes();
    }

    /**
     * @brief Gets the coordinates of the particle's smoothing domain vertices.
     * @return A const reference to an Eigen matrix containing the vertex coordinates.
     */
    const Eigen::Matrix<double, nDim, nVertices>& getSmoothingVertexCoordinates() const
    {
      return _cellForSmoothing.nodes();
    }

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
     * @brief Gets the coordinates of a specific face's center in the intermediate configuration.
     *        This refers to the faces of the *deformed geometry*.
     * @param faceID The ID of the face.
     * @return An Eigen vector representing the face center coordinates.
     */
    Eigen::Matrix<double, nDim, 1> getFaceCenterCoordinates( int faceID ) const
    {
      return _cellForGeometryDeformed.getFaceCenterCoordinates( faceID );
    }

    /**
     * @brief Gets the coordinates of a specific evaluation point (face center of the smoothing domain).
     * @param faceID The ID of the face (1-based).
     * @return An Eigen vector representing the evaluation point coordinate.
     */
    Eigen::Matrix<double, nDim, 1> getSmoothingDomainFaceCenterCoordinates(int faceID) const
    {
        return _cellForSmoothing.getFaceCenterCoordinates(faceID);
    }

    /**
     * @brief Gets the number of evaluation points (number of faces of the smoothing domain).
     * @return The number of evaluation points.
     */
    int getNumberOfFaces() const { return _cellForGeometryDeformed.getNumberOfFaces(); }

    /**
     * @brief Get the smoothing volume of the particle.
     * @return The smoothing volume of the particle.
     */
    double getSmoothingVolume() const { return _cellForSmoothing.volume(); }

    /**
     * @brief Gets the center coordinates of the particle in the intermediate configuration.
     * @return An Eigen vector representing the center coordinates.
     */
    Eigen::Matrix<double, nDim, 1> getCenterCoordinates() const
    {
      return _cellForGeometryDeformed.centroid();
    }

    /**
     * @brief Gets the undeformed volume of the particle.
     * @return The undeformed volume.
     */
    double getVolumeUndeformed() const { return _cellForGeometryUndeformed.volume(); }

    /**
     * @brief Provides a const reference to the vertex displacements of the smoothing domain.
     * @return A const reference to an Eigen matrix containing the vertex displacements.
     */
    const Eigen::Matrix<double, nDim, nVertices>& getSmoothingDomainVertexDisplacements() const
    {
      return _vertex_displacements_smoothingDomain;
    }

    /**
     * @brief Returns the boundary surface vector for a given face ID, from the *deformed geometry*.
     *        This is typically used for distributed loads.
     * @param faceID The ID of the face.
     * @return An Eigen vector representing the boundary surface vector.
     */
    Eigen::Matrix<double, nDim, 1> getFaceBoundaryVector(int faceID) const
    {
        return _cellForGeometryDeformed.boundarySurfaceVector(faceID);
    }

  
  inline std::vector<int> getSubCellIndicesOnParentFace( int parentFaceId ) const
  {
        return _cellForGeometryUndeformed.getSubCellIndicesOnParentFace( parentFaceId );
  }

    /**
     * @brief Returns the boundary surface vector for a given face ID, from the *smoothing domain*.
     *        This is typically used for computing smoothed shape function gradients.
     * @param faceID The ID of the face.
     * @return An Eigen vector representing the boundary surface vector.
     */
    Eigen::Matrix<double, nDim, 1> getSmoothingBoundarySurfaceVector(int faceID) const
    {
        return _cellForSmoothing.boundarySurfaceVector(faceID);
    }

    /**
     * @brief Returns the second moments of area/volume for the deformed geometry.
     * @return An Eigen matrix representing the second moments.
     */
    Eigen::Matrix<double, nDim, nDim> getGeometrySecondMoments() const
    {
        return _cellForGeometryDeformed.secondMoments();
    }

    // Public Methods (Setters/Updaters)
    /**
     * @brief Updates the particle's position and volume to the reference intermediate configuration.
     * @param F_physics The physics-specific deformation gradient.
     * @param centerDisplacement The uniform displacement of the particle's center.
     */
    void acceptStateAndPosition( const Eigen::Matrix< double, nDim, nDim >& F_physics,
                                 const Eigen::Matrix< double, nDim, 1 >&    centerDisplacement)
    {
      _cellForGeometryDeformed.updateVertexCoordinates( _cellForGeometryUndeformed.nodes());
      _cellForGeometryDeformed.applyDeformationGradient( F_physics );
      _cellForGeometryDeformed.applyUniformDisplacement( centerDisplacement );

      const auto FSmoothing = _computeSmoothingDomainDeformationTensorTotal( F_physics );

      _cellForSmoothing.updateVertexCoordinates( _cellForGeometryUndeformed.nodes() );
      _cellForSmoothing.applyDeformationGradient( FSmoothing );
      _cellForSmoothing.applyUniformDisplacement( centerDisplacement );

      _vertex_displacements_smoothingDomain = _cellForSmoothing.nodes() - _cellForGeometryUndeformed.nodes();
    }

    /**
     * @brief Uniformly subdivides the particle domain into smaller domains.
     * @param levels The number of subdivision levels.
     * @return A vector of ParticleDomain instances representing the subdivided particles.
     */
    std::vector< ParticleDomain<nDim, nVertices> > uniformSubdivided() const
    {
        std::vector< ParticleDomain<nDim, nVertices> > subdividedDomains;

        auto subdividedCells = _cellForGeometryUndeformed.uniformSubdivided();

        for (const auto& cell : subdividedCells)
        {
            ParticleDomain<nDim, nVertices> newDomain(
                cell.nodes().data(),
                nDim * nVertices,
                smoothingVolumeUpdateType
            );
            subdividedDomains.push_back(newDomain);
        }

        return subdividedDomains;
    }

    // Public Static Methods
    /**
     * @brief Computes the centroid coordinates from a given set of vertex coordinates.
     * @param vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return An Eigen vector representing the centroid coordinates.
     */
    static CoordinatesSized getCenterFromVertices(
      const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates )
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.centroid();
    }

    /**
     * @brief Computes the volume of a cell defined by a given set of vertex coordinates.
     * @param vertexCoordinates An Eigen matrix containing the vertex coordinates.
     * @return The computed volume.
     */
    static double getVolumeFromVertices( const Eigen::Matrix< double, nDim, nVertices >& vertexCoordinates )
    {
      MarmotLagrangeCell< nDim, nVertices > _lagrangeCell( vertexCoordinates );
      return _lagrangeCell.volume();
    }

  protected:
    // Protected Member Variables
    MarmotLagrangeCell< nDim, nVertices > _cellForGeometryUndeformed; ///< Cell representing the undeformed geometry.
    MarmotLagrangeCell< nDim, nVertices > _cellForGeometryDeformed;   ///< Cell representing the intermediate geometry.
    MarmotLagrangeCell< nDim, nVertices > _cellForSmoothing;          ///< Cell representing the smoothing domain.

    Eigen::Matrix< double, nDim, nVertices > _vertex_displacements_smoothingDomain; ///< Displacements of vertices in the smoothing domain.

    // Protected Methods
    /**
     * @brief Computes the total deformation tensor for the smoothing domain based on the configured update type.
     * @param F_physics The physics-specific deformation gradient.
     * @return An Eigen matrix representing the deformation tensor for the smoothing domain.
     */
    Eigen::Matrix< double, nDim, nDim > _computeSmoothingDomainDeformationTensorTotal(
      const Eigen::Matrix< double, nDim, nDim >& F_physics ) const
    {
      Eigen::Matrix< double, nDim, nDim > F_smoothing;

      switch ( smoothingVolumeUpdateType ) {
      case SmoothingDomainUpdateType::None:
        F_smoothing.setIdentity();
        break;
      case SmoothingDomainUpdateType::DeformationGradient:
        F_smoothing = F_physics;
        break;
      case SmoothingDomainUpdateType::RotationOnly: {
        Eigen::JacobiSVD< Eigen::MatrixXd > svd;
        svd.compute( F_physics, Eigen::ComputeFullU | Eigen::ComputeFullV );
        F_smoothing = svd.matrixU() * svd.matrixV().transpose();
        break;
      }
      case SmoothingDomainUpdateType::RotationAndPrincipalStretch: {
        Eigen::JacobiSVD< Eigen::MatrixXd > svd;
        svd.compute( F_physics, Eigen::ComputeFullU | Eigen::ComputeFullV );
        Eigen::Matrix< double, nDim, nDim > R = svd.matrixU() * svd.matrixV().transpose();
        Eigen::Matrix< double, nDim, nDim > U_ = Eigen::Matrix< double, nDim, nDim >::Identity();
        U_.diagonal()                          = ( R.transpose() * F_physics ).diagonal();
        F_smoothing = R * U_;
        break;
      }
      }
      return F_smoothing;
    }
  };

  template < int nDim, int nVertices >
  ParticleDomain< nDim, nVertices >::ParticleDomain(
    const double*                      vertexCoordinates,
    int                                nVertexCoordinates,
    const SmoothingDomainUpdateType    smoothingVolumeUpdateType )
    : _cellForGeometryUndeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForGeometryDeformed( vertexCoordinates, nVertexCoordinates ),
      _cellForSmoothing( vertexCoordinates, nVertexCoordinates ),
      smoothingVolumeUpdateType( smoothingVolumeUpdateType )
      // _vertexCoordinates_Undeformed( Eigen::Map< const Eigen::Matrix< double, nDim, nVertices > >( vertexCoordinates ) )
  {
  }

} // namespace Marmot::Meshfree
