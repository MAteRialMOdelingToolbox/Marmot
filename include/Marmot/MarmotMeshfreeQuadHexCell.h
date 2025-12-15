#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

/**
 * @file MarmotMeshfreeQuadHexCell.h
 * @brief Defines generic cells (Quad4, Hex8) for meshfree methods.
 *
 * This file provides template classes for 2D quadrilateral (Quad4) and 3D hexahedral (Hex8)
 * isoparametric cells. It includes functionality for shape functions, Jacobian computation,
 * centroid, volume/area, second moments, and uniform subdivision.
 *
 * The elements use Abaqus node and face numbering conventions where applicable.
 */

namespace Marmot::Meshfree {

  /**
   * @brief Generic Gauss rule definition.
   * @tparam nDim The dimension of the space (e.g., 2 for 2D, 3 for 3D).
   *
   * This template struct provides static methods to retrieve Gauss points and weights
   * for numerical integration. Specializations are provided for 2D and 3D.
   */
  template < int nDim >
  struct GaussRule;

  /**
   * @brief Specialization of GaussRule for 2D (2x2 Gauss points).
   */
  template <>
  struct GaussRule< 2 > {
    using Scalar = double; ///< Type for scalar values.
    using Vec    = Eigen::Matrix< Scalar, 2, 1 >; ///< Type for 2D vectors.

    /**
     * @brief Returns the 2x2 Gauss points for a 2D element.
     * @return An array of 4 Eigen::Vector2d representing the Gauss points.
     */
    static std::array< Vec, 4 > points()
    {
      const Scalar a = 1.0 / std::sqrt( 3.0 );
      return { Vec( -a, -a ), Vec( +a, -a ), Vec( +a, +a ), Vec( -a, +a ) };
    }

    /**
     * @brief Returns the weights for the 2x2 Gauss points in 2D.
     * @return An array of 4 doubles, all 1.0 for a 2x2 rule.
     */
    static std::array< Scalar, 4 > weights() { return { 1.0, 1.0, 1.0, 1.0 }; }
  };

  /**
   * @brief Specialization of GaussRule for 3D (2x2x2 Gauss points).
   */
  template <>
  struct GaussRule< 3 > {
    using Scalar = double; ///< Type for scalar values.
    using Vec    = Eigen::Matrix< Scalar, 3, 1 >; ///< Type for 3D vectors.

    /**
     * @brief Returns the 2x2x2 Gauss points for a 3D element.
     * @return An array of 8 Eigen::Vector3d representing the Gauss points.
     */
    static std::array< Vec, 8 > points()
    {
      const Scalar         a = 1.0 / std::sqrt( 3.0 );
      std::array< Vec, 8 > gps;
      int                  k = 0;
      for ( Scalar xi : { -a, +a } )
        for ( Scalar eta : { -a, +a } )
          for ( Scalar z : { -a, +a } ) {
            Vec s;
            s << xi, eta, z;
            gps[k++] = s;
          }
      return gps;
    }

    /**
     * @brief Returns the weights for the 2x2x2 Gauss points in 3D.
     * @return An array of 8 doubles, all 1.0 for a 2x2x2 rule.
     */
    static std::array< Scalar, 8 > weights() { return { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }; }
  };

  /**
   * @brief Generic isoparametric Lagrange element.
   * @tparam nDim The dimension of the element (e.g., 2 for Quad4, 3 for Hex8).
   * @tparam nNodes The number of nodes in the element (e.g., 4 for Quad4, 8 for Hex8).
   *
   * This class provides common functionality for isoparametric elements, including
   * mapping from natural to physical coordinates, Jacobian computation, and
   * integration-based properties like centroid and volume.
   *
   * Specific shape functions and derivatives must be specialized for concrete topologies.
   */
  template < int nDim, int nNodes >
  class MarmotLagrangeCell {
  public:
    using Scalar = double; ///< Type for scalar values.
    using Vec    = Eigen::Matrix< Scalar, nDim, 1 >; ///< Type for vectors in nDim.
    using Mat    = Eigen::Matrix< Scalar, nDim, nNodes >; ///< Type for node coordinate matrix.
    using RowN   = Eigen::Matrix< Scalar, 1, nNodes >; ///< Type for shape function vector.
    using dNMat  = Eigen::Matrix< Scalar, nDim, nNodes >; ///< Type for shape function derivative matrix.

    using Self = MarmotLagrangeCell< nDim, nNodes >; ///< Alias for the current class type.

    /**
     * @brief Default constructor.
     */
    MarmotLagrangeCell() = default;

    /**
     * @brief Constructor that initializes the element with node coordinates.
     * @param nodes An Eigen matrix where each column represents a node's coordinates.
     */
    explicit MarmotLagrangeCell( const Mat& nodes ) : _nodes( nodes ) {}

    /**
     * @brief Constructor that initializes the element with raw node data.
     * @param nodesData A pointer to a flat array of node coordinates (e.g., [x1,y1,z1, x2,y2,z2, ...]).
     * @param nNodesData The total number of scalar values in nodesData (nDim * nNodes).
     * @throws std::invalid_argument if nNodesData does not match nDim * nNodes.
     */
    explicit MarmotLagrangeCell( const Scalar* nodesData, int nNodesData )
    {
      if ( nNodesData != nDim * nNodes )
        throw std::invalid_argument( "MarmotLagrangeCell: invalid nodes data size" );
      _nodes = Eigen::Map< const Mat >( nodesData, nDim, nNodes );

      // std::cout << _nodes << std::endl;
      // std::cout << "MarmotLagrangeCell created with " << nNodes << " nodes in " << nDim
      //           << "D." << std::endl;
    }

    /**
     * @brief Returns a const reference to the node coordinate matrix.
     * @return A const Eigen matrix of node coordinates.
     */
    const Mat& nodes() const { return _nodes; }

    /**
     * @brief Returns a mutable reference to the node coordinate matrix.
     * @return A mutable Eigen matrix of node coordinates.
     */
    Mat& nodes() { return _nodes; }

    // ------------------------------------------------------------------------
    // Element-specific: shape functions and derivatives (must be specialized)
    // ------------------------------------------------------------------------
    /**
     * @brief Computes the shape functions at a given natural coordinate.
     * @param s The natural coordinate vector (xi, eta, zeta).
     * @return A row vector of shape function values.
     * @note This method must be specialized for each concrete element type.
     */
    static RowN shapeFunctions( const Vec& s );

    /**
     * @brief Computes the derivatives of the shape functions with respect to natural coordinates.
     * @param s The natural coordinate vector (xi, eta, zeta).
     * @return A matrix where rows are derivatives w.r.t. natural coordinates and columns are nodes.
     * @note This method must be specialized for each concrete element type.
     */
    static dNMat dShapeFunctions( const Vec& s );

    /**
     * @brief Returns the interpolation operator (same as shape functions for Lagrange elements).
     * @param s The natural coordinate vector.
     * @return A row vector of shape function values.
     */
    static RowN interpolationOperator( const Vec& s ) { return shapeFunctions( s ); }

    // ------------------------------------------------------------------------
    // Mapping and Jacobian
    // ------------------------------------------------------------------------
    /**
     * @brief Maps a natural coordinate to a physical coordinate.
     * @param s The natural coordinate vector.
     * @return The physical coordinate vector.
     */
    Vec mapToPhysical( const Vec& s ) const { return _nodes * shapeFunctions( s ).transpose(); }

    /**
     * @brief Computes the Jacobian matrix of the mapping from natural to physical coordinates.
     * @param s The natural coordinate vector.
     * @return The Jacobian matrix (nDim x nDim).
     */
    Eigen::Matrix< Scalar, nDim, nDim > jacobian( const Vec& s ) const
    {
      dNMat dN = dShapeFunctions( s );
      return _nodes * dN.transpose();
    }

    // ------------------------------------------------------------------------
    // Generic exact centroid via Gauss integration
    // ------------------------------------------------------------------------
    /**
     * @brief Computes the centroid of the element using Gauss integration.
     * @return The physical coordinate vector of the centroid.
     */
    Vec centroid() const
    {
      Vec    c   = Vec::Zero();
      Scalar V   = 0.0;
      auto   gps = GaussRule< nDim >::points();
      auto   wgs = GaussRule< nDim >::weights();

      for ( size_t i = 0; i < gps.size(); ++i ) {
        const auto&                         s    = gps[i];
        const auto&                         w    = wgs[i];
        dNMat                               dN   = dShapeFunctions( s );
        Eigen::Matrix< Scalar, nDim, nDim > J    = _nodes * dN.transpose();
        const Scalar                        detJ = std::abs( J.determinant() ) * w;

        Vec x = mapToPhysical( s );

        V += detJ;
        c += x * detJ;
      }

      return c / V;
    }

    // ------------------------------------------------------------------------
    // Generic exact volume/area via Gauss integration
    // ------------------------------------------------------------------------
    /**
     * @brief Computes the volume (or area in 2D) of the element using Gauss integration.
     * @return The scalar volume/area.
     */
    Scalar volume() const
    {
      Scalar V   = 0.0;
      auto   gps = GaussRule< nDim >::points();
      auto   wgs = GaussRule< nDim >::weights();

      for ( size_t i = 0; i < gps.size(); ++i ) {
        const auto& s = gps[i];
        const auto& w = wgs[i];

        dNMat                               dN   = dShapeFunctions( s );
        Eigen::Matrix< Scalar, nDim, nDim > J    = _nodes * dN.transpose();
        const Scalar                        detJ = std::abs( J.determinant() ) * w;
        V += detJ;
      }
      return V;
    }

    /**
     * @brief Computes the area of the element (only valid for 2D elements).
     * @return The scalar area.
     * @note This method asserts that nDim is 2.
     */
    Scalar area() const
    {
      static_assert( nDim == 2, "area() only makes sense in 2D." );
      return volume();
    }

    // ------------------------------------------------------------------------
    // Full second moment matrix: I_kj = ∫ d_k d_j dV
    // d = x - C
    // ------------------------------------------------------------------------
    /**
     * @brief Computes the second moment matrix of the element about its centroid.
     *
     * The second moment matrix I_kj is defined as ∫ (x_k - C_k)(x_j - C_j) dV,
     * where C is the centroid and dV is the differential volume.
     *
     * @return An Eigen matrix representing the second moment matrix.
     */
    Eigen::Matrix< Scalar, nDim, nDim > secondMoments() const
    {
      Vec C = centroid();

      Eigen::Matrix< Scalar, nDim, nDim > I;
      I.setZero();

      auto gps = GaussRule< nDim >::points();
      auto wgs = GaussRule< nDim >::weights();

      for ( size_t i = 0; i < gps.size(); ++i ) {
        const auto& s = gps[i];
        const auto& w = wgs[i];

        dNMat                               dN   = dShapeFunctions( s );
        Eigen::Matrix< Scalar, nDim, nDim > J    = _nodes * dN.transpose();
        Scalar                              detJ = std::abs( J.determinant() ) * w;

        // physical coordinates
        Vec x = mapToPhysical( s );
        Vec d = x - C; // deviation from centroid

        // Accumulate outer product: d_k d_j
        I += ( d * d.transpose() ) * detJ;
      }

      return I;
    }

    // ------------------------------------------------------------------------
    // Update vertex coordinates
    // ------------------------------------------------------------------------
    /**
     * @brief Updates the coordinates of the element's nodes.
     * @param newNodes The new Eigen matrix of node coordinates.
     */
    void updateVertexCoordinates( const Mat& newNodes ) { _nodes = newNodes; }

    // ------------------------------------------------------------------------
    // Apply deformation gradient relative to centroid:
    // x_new = C + F * (x_old - C)
    // ------------------------------------------------------------------------
    /**
     * @brief Applies a deformation gradient relative to the element's centroid.
     *
     * Each node's new position is calculated as: x_new = C + F * (x_old - C),
     * where C is the current centroid and F is the deformation gradient.
     *
     * @param F The deformation gradient matrix.
     */
    void applyDeformationGradient( const Eigen::Matrix< Scalar, nDim, nDim >& F )
    {
      Vec C = centroid();
      for ( int i = 0; i < nNodes; ++i ) {
        Vec Xrel        = _nodes.col( i ) - C;
        _nodes.col( i ) = C + F * Xrel;
      }
    }

    /**
     * @brief Applies a uniform displacement to all nodes of the element.
     * @param dX The displacement vector.
     */
    void applyUniformDisplacement( const Vec& dX )
    {
      for ( int i = 0; i < nNodes; ++i ) {
        _nodes.col( i ) += dX;
      }
    }

    // ------------------------------------------------------------------------
    // Faces: center coords and boundary surface vectors
    //   - must be specialized for concrete topologies
    // ------------------------------------------------------------------------
    /**
     * @brief Computes the center coordinates of a specified face.
     * @param faceId The ID of the face (1-based, typically Abaqus convention).
     * @return The physical coordinate vector of the face center.
     * @throws std::logic_error if not implemented for the specific element type.
     * @note This method must be specialized for each concrete element type.
     */
    Vec getFaceCenterCoordinates( int faceId ) const
    {
      throw std::logic_error( "getFaceCenterCoordinates not implemented for this (nDim,nNodes)" );
    }

    /**
     * @brief Computes the boundary surface vector (normal scaled by area/length) for a face.
     * @param faceId The ID of the face (1-based, typically Abaqus convention).
     * @return The boundary surface vector.
     * @throws std::logic_error if not implemented for the specific element type.
     * @note This method must be specialized for each concrete element type.
     */
    Vec boundarySurfaceVector( int faceId ) const
    {
      throw std::logic_error( "boundarySurfaceVector not implemented for this (nDim,nNodes)" );
    }

    // ------------------------------------------------------------------------
    // Uniform subdivision: must be specialized per topology
    // ------------------------------------------------------------------------
    /**
     * @brief Subdivides the element uniformly into smaller elements.
     * @return A vector of new MarmotLagrangeCell instances representing the subdivided elements.
     * @throws std::logic_error if not implemented for the specific element type.
     * @note This method must be specialized for each concrete element type.
     */
    std::vector< Self > uniformSubdivided( ) const
    {
      // The actual subdivision logic is handled by the specialized versions below.
      throw std::logic_error( "uniformSubdivided not implemented for this (nDim,nNodes)" );
    }

    /**
     * @brief Retrieves the indices of sub-cells that lie on a specified face of the parent cell
     *        after uniform subdivision.
     * @param parentFaceId The ID of the parent cell's face (1-based, Abaqus convention).
     * @return A vector of indices of the sub-cells in the flat list returned by `uniformSubdivided`.
     * @throws std::out_of_range if parentFaceId is invalid.
     * @throws std::logic_error if not implemented for the specific element type.
     */
    std::vector<int> getSubCellIndicesOnParentFace( int parentFaceId) const;

    // ------------------------------------------------------------------------
    // Ensight Gold cell shape name (specialized for supported element types)
    // ------------------------------------------------------------------------
    /**
     * @brief Returns the Ensight Gold cell shape name for the element.
     * @return A string representing the Ensight Gold cell shape (e.g., "quad4", "hexa8").
     * @throws std::logic_error if not implemented for the specific element type.
     * @note This method must be specialized for each concrete element type.
     */
    std::string getCellShape() const
    {
      throw std::logic_error( "getCellShape() not implemented for this (nDim,nNodes)" );
    }

    /**
     * @brief Returns the number of nodes in the element.
     * @return The integer number of nodes.
     */
    int getNumberOfNodes() const { return nNodes; }

    /**
     * @brief Returns the number of dimensions of the element.
     * @return The integer number of dimensions.
     */
    int getNumberOfDimensions() const { return nDim; }

    /**
     * @brief Returns the number of faces of the element.
     * @return The integer number of faces.
     * @throws std::logic_error if not implemented for the specific element type.
     * @note This method must be specialized for each concrete element type.
     */
    int getNumberOfFaces() const
    {
      throw std::logic_error( "getNumberOfFaces() not implemented for this (nDim,nNodes)" );
    }

  private:
    Mat _nodes; ///< Matrix storing the physical coordinates of the element's nodes.

    // Helper to get global grid coordinates for a sub-cell
    // This function determines the (i,j) or (i,j,k) grid position
    // of a sub-cell given its flat index in the `uniformSubdivided` vector
    // The coordinates are 0-indexed, from 0 to (2^levels - 1) along each dimension.
    // NOTE: This helper is now only relevant for levels=1, simplifying its usage.
    std::array<int, nDim> getGlobalGridCoords(int cell_idx ) const
    {

      std::array<int, nDim> coords;
      coords.fill( 0 );

      // The cell_idx directly maps to the local (i,j) or (i,j,k) index.
      int local_child_idx         = cell_idx; // For levels=1, cell_idx is the direct index

      if ( nDim == 2 ) {
        // Quad4 specific local mapping: (0,0), (1,0), (1,1), (0,1)
        int local_i, local_j;
        if ( local_child_idx == 0 ) { // Bottom-left
          local_i = 0;
          local_j = 0;
        }
        else if ( local_child_idx == 1 ) { // Bottom-right
          local_i = 1;
          local_j = 0;
        }
        else if ( local_child_idx == 2 ) { // Top-right
          local_i = 1;
          local_j = 1;
        }
        else { // local_child_idx == 3 (Top-left)
          local_i = 0;
          local_j = 1;
        }
        coords[0] = local_i;
        coords[1] = local_j;
      }
      else { // nDim == 3
        // Hex8 specific local mapping: (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1)
        int local_i, local_j, local_k;
        if ( local_child_idx == 0 ) { // bottom-front-left
          local_i = 0; local_j = 0; local_k = 0;
        }
        else if ( local_child_idx == 1 ) { // bottom-front-right
          local_i = 1; local_j = 0; local_k = 0;
        }
        else if ( local_child_idx == 2 ) { // bottom-back-right
          local_i = 1; local_j = 1; local_k = 0;
        }
        else if ( local_child_idx == 3 ) { // bottom-back-left
          local_i = 0; local_j = 1; local_k = 0;
        }
        else if ( local_child_idx == 4 ) { // top-front-left
          local_i = 0; local_j = 0; local_k = 1;
        }
        else if ( local_child_idx == 5 ) { // top-front-right
          local_i = 1; local_j = 0; local_k = 1;
        }
        else if ( local_child_idx == 6 ) { // top-back-right
          local_i = 1; local_j = 1; local_k = 1;
        }
        else { // local_child_idx == 7 (top-back-left)
          local_i = 0; local_j = 1; local_k = 1;
        }
        coords[0] = local_i;
        coords[1] = local_j;
        coords[2] = local_k;
      }
      return coords;
    }
  };

  // ============================================================================
  //   Specialization: Quad4  (nDim=2, nNodes=4) – Abaqus node numbering
  // ============================================================================

  /**
   * @brief Specialization of shapeFunctions for a 2D, 4-node quadrilateral (Quad4).
   * @param s The natural coordinate vector (xi, eta).
   * @return A row vector of shape function values for Quad4.
   */
  template <>
  inline MarmotLagrangeCell< 2, 4 >::RowN MarmotLagrangeCell< 2, 4 >::shapeFunctions( const Vec& s )
  {
    Scalar xi = s[0], eta = s[1];
    RowN   N;
    N( 0 ) = 0.25 * ( 1 - xi ) * ( 1 - eta ); // node 1
    N( 1 ) = 0.25 * ( 1 + xi ) * ( 1 - eta ); // node 2
    N( 2 ) = 0.25 * ( 1 + xi ) * ( 1 + eta ); // node 3
    N( 3 ) = 0.25 * ( 1 - xi ) * ( 1 + eta ); // node 4
    return N;
  }

  /**
   * @brief Specialization of dShapeFunctions for a 2D, 4-node quadrilateral (Quad4).
   * @param s The natural coordinate vector (xi, eta).
   * @return A matrix of shape function derivatives for Quad4.
   */
  template <>
  inline MarmotLagrangeCell< 2, 4 >::dNMat MarmotLagrangeCell< 2, 4 >::dShapeFunctions( const Vec& s )
  {
    Scalar xi = s[0], eta = s[1];
    dNMat  dN;

    // d/dxi
    dN( 0, 0 ) = -0.25 * ( 1 - eta );
    dN( 0, 1 ) = 0.25 * ( 1 - eta );
    dN( 0, 2 ) = 0.25 * ( 1 + eta );
    dN( 0, 3 ) = -0.25 * ( 1 + eta );

    // d/deta
    dN( 1, 0 ) = -0.25 * ( 1 - xi );
    dN( 1, 1 ) = -0.25 * ( 1 + xi );
    dN( 1, 2 ) = 0.25 * ( 1 + xi );
    dN( 1, 3 ) = 0.25 * ( 1 - xi );

    return dN;
  }

  /**
   * @brief Specialization of getFaceCenterCoordinates for Quad4.
   * @param faceId The ID of the face (1-4, Abaqus convention).
   * @return The physical coordinate vector of the face center.
   * @throws std::out_of_range if faceId is not between 1 and 4.
   */
  template <>
  inline MarmotLagrangeCell< 2, 4 >::Vec MarmotLagrangeCell< 2, 4 >::getFaceCenterCoordinates( int faceId ) const
  {
    static const int edges[4][2] = {
      { 0, 1 }, // S1 (Abaqus)
      { 1, 2 }, // S2 (Abaqus)
      { 2, 3 }, // S3 (Abaqus)
      { 3, 0 }  // S4 (Abaqus)
    };

    if ( faceId < 1 || faceId > 4 )
      throw std::out_of_range( "Quad4 faceId must be 1..4 (Abaqus)" );

    int i0 = edges[faceId - 1][0];
    int i1 = edges[faceId - 1][1];
    return 0.5 * ( _nodes.col( i0 ) + _nodes.col( i1 ) );
  }

  /**
   * @brief Specialization of boundarySurfaceVector for Quad4.
   *
   * Computes the outward normal vector scaled by the edge length for a given face.
   * Assumes counter-clockwise node ordering for outward normal.
   *
   * @param faceId The ID of the face (1-4, Abaqus convention).
   * @return The boundary surface vector.
   * @throws std::out_of_range if faceId is not between 1 and 4.
   */
  template <>
  inline MarmotLagrangeCell< 2, 4 >::Vec MarmotLagrangeCell< 2, 4 >::boundarySurfaceVector( int faceId ) const
  {
    static const int edges[4][2] = {
      { 0, 1 }, // S1 (Abaqus)
      { 1, 2 }, // S2 (Abaqus)
      { 2, 3 }, // S3 (Abaqus)
      { 3, 0 }  // S4 (Abaqus)
    };

    if ( faceId < 1 || faceId > 4 )
      throw std::out_of_range( "Quad4 faceId must be 1..4 (Abaqus)" );

    int i0 = edges[faceId - 1][0];
    int i1 = edges[faceId - 1][1];

    Vec e = _nodes.col( i1 ) - _nodes.col( i0 );
    return Vec( e[1], -e[0] ); // outward normal * edge length (for CCW nodes)
  }

  /**
   * @brief Specialization of uniformSubdivided for Quad4.
   * @return A vector of new Quad4 elements.
   */
  template <>
  inline std::vector< MarmotLagrangeCell< 2, 4 > > MarmotLagrangeCell< 2, 4 >::uniformSubdivided( ) const
  {

    std::vector< Self > next;
    next.reserve( 4 );

    Vec s;
    Vec grid[3][3]; // Stores physical coordinates of points at -1, 0, +1 in natural space

    for ( int j = 0; j < 3; ++j)
      for ( int i = 0; i < 3; ++i ) {
        s << -1.0 + i, -1.0 + j; // Natural coordinates (-1,-1), (0,-1), (1,-1), etc.
        grid[i][j] = mapToPhysical( s );
      }

    auto makeChild = [&]( int i0, int j0 ) {
      Mat m;
      m.col( 0 ) = grid[i0][j0];
      m.col( 1 ) = grid[i0 + 1][j0];
      m.col( 2 ) = grid[i0 + 1][j0 + 1];
      m.col( 3 ) = grid[i0][j0 + 1];
      return Self( m );
    };

    next.push_back( makeChild( 0, 0 ) ); // bottom-left
    next.push_back( makeChild( 1, 0 ) ); // bottom-right
    next.push_back( makeChild( 1, 1 ) ); // top-right
    next.push_back( makeChild( 0, 1 ) ); // top-left

    return next;
  }

  /**
   * @brief Specialization of getSubCellIndicesOnParentFace for Quad4.
   * @param parentFaceId The ID of the parent cell's face (1-4, Abaqus convention).
   * @return A vector of indices of the sub-cells in the flat list returned by `uniformSubdivided`.
   * @throws std::out_of_range if parentFaceId is not between 1 and 4.
   */
  template <>
  inline std::vector<int> MarmotLagrangeCell< 2, 4 >::getSubCellIndicesOnParentFace( int parentFaceId ) const
  {

    std::vector<int> indices;

    switch ( parentFaceId ) {
        case 1: // S1: Bottom (eta=-1) -> sub-cells 0 and 1 are on this face
            indices.push_back(0);
            indices.push_back(1);
            break;
        case 2: // S2: Right (xi=+1) -> sub-cells 1 and 2 are on this face
            indices.push_back(1);
            indices.push_back(2);
            break;
        case 3: // S3: Top (eta=+1) -> sub-cells 2 and 3 are on this face
            indices.push_back(2);
            indices.push_back(3);
            break;
        case 4: // S4: Left (xi=-1) -> sub-cells 0 and 3 are on this face
            indices.push_back(0);
            indices.push_back(3);
            break;
        default:
            throw std::out_of_range( "Quad4 parentFaceId must be 1..4 (Abaqus)" );
    }
    return indices;
  }

  /**
   * @brief Specialization of getCellShape for Quad4.
   * @return The string "quad4" (Ensight Gold notation).
   */
  template <>
  inline std::string MarmotLagrangeCell< 2, 4 >::getCellShape() const
  {
    return "quad4"; // Ensight Gold notation
  }

  /**
   * @brief Specialization of getNumberOfFaces for Quad4.
   * @return The integer 4.
   */
  template <>
  inline int MarmotLagrangeCell< 2, 4 >::getNumberOfFaces() const
  {
    return 4;
  }

  // ============================================================================
  //   Specialization: Hex8  (nDim=3, nNodes=8) – Abaqus node & face numbering
  // ============================================================================

  /**
   * @brief Specialization of shapeFunctions for a 3D, 8-node hexahedron (Hex8).
   * @param s The natural coordinate vector (xi, eta, zeta).
   * @return A row vector of shape function values for Hex8.
   */
  template <>
  inline MarmotLagrangeCell< 3, 8 >::RowN MarmotLagrangeCell< 3, 8 >::shapeFunctions( const Vec& s )
  {
    Scalar xi = s[0], eta = s[1], z = s[2];
    Scalar c = 0.125;

    RowN N;
    N( 0 ) = c * ( 1 - xi ) * ( 1 - eta ) * ( 1 - z );
    N( 1 ) = c * ( 1 + xi ) * ( 1 - eta ) * ( 1 - z );
    N( 2 ) = c * ( 1 + xi ) * ( 1 + eta ) * ( 1 - z );
    N( 3 ) = c * ( 1 - xi ) * ( 1 + eta ) * ( 1 - z );
    N( 4 ) = c * ( 1 - xi ) * ( 1 - eta ) * ( 1 + z );
    N( 5 ) = c * ( 1 + xi ) * ( 1 - eta ) * ( 1 + z );
    N( 6 ) = c * ( 1 + xi ) * ( 1 + eta ) * ( 1 + z );
    N( 7 ) = c * ( 1 - xi ) * ( 1 + eta ) * ( 1 + z );
    return N;
  }

  /**
   * @brief Specialization of dShapeFunctions for a 3D, 8-node hexahedron (Hex8).
   * @param s The natural coordinate vector (xi, eta, zeta).
   * @return A matrix of shape function derivatives for Hex8.
   */
  template <>
  inline MarmotLagrangeCell< 3, 8 >::dNMat MarmotLagrangeCell< 3, 8 >::dShapeFunctions( const Vec& s )
  {
    Scalar xi = s[0], eta = s[1], z = s[2];
    Scalar c = 0.125;
    dNMat  dN;

    // d/dxi
    dN( 0, 0 ) = -c * ( 1 - eta ) * ( 1 - z );
    dN( 0, 1 ) = c * ( 1 - eta ) * ( 1 - z );
    dN( 0, 2 ) = c * ( 1 + eta ) * ( 1 - z );
    dN( 0, 3 ) = -c * ( 1 + eta ) * ( 1 - z );
    dN( 0, 4 ) = -c * ( 1 - eta ) * ( 1 + z );
    dN( 0, 5 ) = c * ( 1 - eta ) * ( 1 + z );
    dN( 0, 6 ) = c * ( 1 + eta ) * ( 1 + z );
    dN( 0, 7 ) = -c * ( 1 + eta ) * ( 1 + z );

    // d/deta
    dN( 1, 0 ) = -c * ( 1 - xi ) * ( 1 - z );
    dN( 1, 1 ) = -c * ( 1 + xi ) * ( 1 - z );
    dN( 1, 2 ) = c * ( 1 + xi ) * ( 1 - z );
    dN( 1, 3 ) = c * ( 1 - xi ) * ( 1 - z );
    dN( 1, 4 ) = -c * ( 1 - xi ) * ( 1 + z );
    dN( 1, 5 ) = -c * ( 1 + xi ) * ( 1 + z );
    dN( 1, 6 ) = c * ( 1 + xi ) * ( 1 + z );
    dN( 1, 7 ) = c * ( 1 - xi ) * ( 1 + z );

    // d/dz
    dN( 2, 0 ) = -c * ( 1 - xi ) * ( 1 - eta );
    dN( 2, 1 ) = -c * ( 1 + xi ) * ( 1 - eta );
    dN( 2, 2 ) = -c * ( 1 + xi ) * ( 1 + eta );
    dN( 2, 3 ) = -c * ( 1 - xi ) * ( 1 + eta );
    dN( 2, 4 ) = c * ( 1 - xi ) * ( 1 - eta );
    dN( 2, 5 ) = c * ( 1 + xi ) * ( 1 - eta );
    dN( 2, 6 ) = c * ( 1 + xi ) * ( 1 + eta );
    dN( 2, 7 ) = c * ( 1 - xi ) * ( 1 + eta );

    return dN;
  }

  /**
   * @brief Specialization of getFaceCenterCoordinates for Hex8.
   * @param faceId The ID of the face (1-6, Abaqus convention).
   * @return The physical coordinate vector of the face center.
   * @throws std::out_of_range if faceId is not between 1 and 6.
   */
  template <>
  inline MarmotLagrangeCell< 3, 8 >::Vec MarmotLagrangeCell< 3, 8 >::getFaceCenterCoordinates( int faceId ) const
  {
    // Abaqus C3D8 Standard Faces (0-based indices)
    static const int faces[6][4] = {
      { 0, 3, 2, 1 }, // S1: Bottom (points -Z local)
      { 4, 5, 6, 7 }, // S2: Top    (points +Z local)
      { 0, 1, 5, 4 }, // S3: Front  (points -Y local)
      { 1, 2, 6, 5 }, // S4: Right  (points +X local)
      { 2, 3, 7, 6 }, // S5: Back   (points +Y local)
      { 3, 0, 4, 7 }  // S6: Left   (points -X local)
    };
    if ( faceId < 1 || faceId > 6 )
      throw std::out_of_range( "Hex8 faceId must be 1..6 (Abaqus)" );

    Vec c = Vec::Zero();
    for ( int i = 0; i < 4; ++i )
      c += _nodes.col( faces[faceId - 1][i] );

    return c / 4.0;
  }

  /**
   * @brief Specialization of boundarySurfaceVector for Hex8.
   *
   * Computes the outward normal vector scaled by the face area for a given face.
   * Uses the cross product of the face diagonals for a robust normal estimation
   * on potentially warped hexahedral faces.
   *
   * @param faceId The ID of the face (1-6, Abaqus convention).
   * @return The boundary surface vector.
   * @throws std::out_of_range if faceId is not between 1 and 6.
   */
  template <>
  inline MarmotLagrangeCell< 3, 8 >::Vec MarmotLagrangeCell< 3, 8 >::boundarySurfaceVector( int faceId ) const
  {
    // Abaqus C3D8 Standard Faces (0-based indices)
    // Note: To get OUTWARD normals, we use the Right Hand Rule counter-clockwise
    // ordering when looking at the face from the OUTSIDE.
    static const int faces[6][4] = {
      { 0, 3, 2, 1 }, // S1: Bottom (points -Z local)
      { 4, 5, 6, 7 }, // S2: Top    (points +Z local)
      { 0, 1, 5, 4 }, // S3: Front  (points -Y local)
      { 1, 2, 6, 5 }, // S4: Right  (points +X local)
      { 2, 3, 7, 6 }, // S5: Back   (points +Y local)
      { 3, 0, 4, 7 }  // S6: Left   (points -X local)
    };

    if ( faceId < 1 || faceId > 6 )
      throw std::out_of_range( "Hex8 faceId must be 1..6 (Abaqus)" );

    const int* f = faces[faceId - 1];

    Vec p1 = _nodes.col( f[0] );
    Vec p2 = _nodes.col( f[1] );
    Vec p3 = _nodes.col( f[2] );
    Vec p4 = _nodes.col( f[3] );

    // Area Vector = 0.5 * (Diagonal A x Diagonal B)
    // Diagonal A = p3 - p1
    // Diagonal B = p4 - p2
    Vec diag1 = p3 - p1;
    Vec diag2 = p4 - p2;

    return 0.5 * diag1.cross( diag2 );
  }

  /**
   * @brief Specialization of uniformSubdivided for Hex8.
   * @return A vector of new Hex8 elements.
   */
  template <>
  inline std::vector< MarmotLagrangeCell< 3, 8 > > MarmotLagrangeCell< 3, 8 >::uniformSubdivided( ) const
  {

    std::vector< Self > next;
    next.reserve( 8 );

    Vec s;
    Vec grid[3][3][3]; // Stores physical coordinates of points at -1, 0, +1 in natural space

    for ( int k = 0; k < 3; ++k )
      for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i ) {
          s << -1.0 + i, -1.0 + j, -1.0 + k; // Natural coordinates (-1,-1,-1), (0,-1,-1), etc.
          grid[i][j][k] = mapToPhysical( s );
        }

    auto makeChild = [&]( int i0, int j0, int k0 ) {
      Mat m;
      m.col( 0 ) = grid[i0][j0][k0];
      m.col( 1 ) = grid[i0 + 1][j0][k0];
      m.col( 2 ) = grid[i0 + 1][j0 + 1][k0];
      m.col( 3 ) = grid[i0][j0 + 1][k0];
      m.col( 4 ) = grid[i0][j0][k0 + 1];
      m.col( 5 ) = grid[i0 + 1][j0][k0 + 1];
      m.col( 6 ) = grid[i0 + 1][j0 + 1][k0 + 1];
      m.col( 7 ) = grid[i0][j0 + 1][k0 + 1];
      return Self( m );
    };

    next.push_back( makeChild( 0, 0, 0 ) ); // bottom-front-left
    next.push_back( makeChild( 1, 0, 0 ) ); // bottom-front-right
    next.push_back( makeChild( 1, 1, 0 ) ); // bottom-back-right
    next.push_back( makeChild( 0, 1, 0 ) ); // bottom-back-left
    next.push_back( makeChild( 0, 0, 1 ) ); // top-front-left
    next.push_back( makeChild( 1, 0, 1 ) ); // top-front-right
    next.push_back( makeChild( 1, 1, 1 ) ); // top-back-right
    next.push_back( makeChild( 0, 1, 1 ) ); // top-back-left

    return next;
  }

  /**
   * @brief Specialization of getSubCellIndicesOnParentFace for Hex8.
   * @param parentFaceId The ID of the parent cell's face (1-6, Abaqus convention).
   * @return A vector of indices of the sub-cells in the flat list returned by `uniformSubdivided`.
   * @throws std::out_of_range if parentFaceId is not between 1 and 6.
   */
  template <>
  inline std::vector<int> MarmotLagrangeCell< 3, 8 >::getSubCellIndicesOnParentFace( int parentFaceId ) const
  {

    std::vector<int> indices;

    switch ( parentFaceId ) {
        case 1: // S1: Bottom (z=-1) -> sub-cells 0, 1, 2, 3 are on this face
            indices.push_back(0);
            indices.push_back(1);
            indices.push_back(2);
            indices.push_back(3);
            break;
        case 2: // S2: Top (z=+1) -> sub-cells 4, 5, 6, 7 are on this face
            indices.push_back(4);
            indices.push_back(5);
            indices.push_back(6);
            indices.push_back(7);
            break;
        case 3: // S3: Front (eta=-1) -> sub-cells 0, 1, 4, 5 are on this face
            indices.push_back(0);
            indices.push_back(1);
            indices.push_back(4);
            indices.push_back(5);
            break;
        case 4: // S4: Right (xi=+1) -> sub-cells 1, 2, 5, 6 are on this face
            indices.push_back(1);
            indices.push_back(2);
            indices.push_back(5);
            indices.push_back(6);
            break;
        case 5: // S5: Back (eta=+1) -> sub-cells 2, 3, 6, 7 are on this face
            indices.push_back(2);
            indices.push_back(3);
            indices.push_back(6);
            indices.push_back(7);
            break;
        case 6: // S6: Left (xi=-1) -> sub-cells 0, 3, 4, 7 are on this face
            indices.push_back(0);
            indices.push_back(3);
            indices.push_back(4);
            indices.push_back(7);
            break;
        default:
            throw std::out_of_range( "Hex8 parentFaceId must be 1..6 (Abaqus)" );
    }
    return indices;
  }

  /**
   * @brief Specialization of getCellShape for Hex8.
   * @return The string "hexa8" (Ensight Gold notation).
   */
  template <>
  inline std::string MarmotLagrangeCell< 3, 8 >::getCellShape() const
  {
    return "hexa8"; // Ensight Gold notation
  }

  /**
   * @brief Specialization of getNumberOfFaces for Hex8.
   * @return The integer 6.
   */
  template <>
  inline int MarmotLagrangeCell< 3, 8 >::getNumberOfFaces() const
  {
    return 6;
  }

  // ============================================================================
  // Convenient aliases
  // ============================================================================

  /**
   * @brief Alias for a 2D, 4-node quadrilateral element.
   */
  using Quad4 = MarmotLagrangeCell< 2, 4 >;

  /**
   * @brief Alias for a 3D, 8-node hexahedral element.
   */
  using Hex8  = MarmotLagrangeCell< 3, 8 >;

} // namespace Marmot::Meshfree
