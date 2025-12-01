#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Marmot::Meshfree {

  // ============================================================================
  //   Gauss rules (generic by dimension, specialized for 2D / 3D)
  //   - 2D: 2x2 Gauss
  //   - 3D: 2x2x2 Gauss
  // ============================================================================

  template < int nDim >
  struct GaussRule;

  template <>
  struct GaussRule< 2 > {
    using Scalar = double;
    using Vec    = Eigen::Matrix< Scalar, 2, 1 >;

    static std::array< Vec, 4 > points()
    {
      const Scalar a = 1.0 / std::sqrt( 3.0 );
      return { Vec( -a, -a ), Vec( +a, -a ), Vec( +a, +a ), Vec( -a, +a ) };
    }
    static std::array< Scalar, 4 > weights() { return { 1.0, 1.0, 1.0, 1.0 }; }
  };

  template <>
  struct GaussRule< 3 > {
    using Scalar = double;
    using Vec    = Eigen::Matrix< Scalar, 3, 1 >;

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

    static std::array< Scalar, 8 > weights() { return { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }; }
  };

  // ============================================================================
  //   Generic isoparametric element: nDim = 2 or 3, nNodes = 4 or 8 here
  //   Specialized for Quad4 (2,4) and Hex8 (3,8)
  // ============================================================================

  template < int nDim, int nNodes >
  class MarmotLagrangeCell {
  public:
    using Scalar = double;
    using Vec    = Eigen::Matrix< Scalar, nDim, 1 >;
    using Mat    = Eigen::Matrix< Scalar, nDim, nNodes >;
    using RowN   = Eigen::Matrix< Scalar, 1, nNodes >;
    using dNMat  = Eigen::Matrix< Scalar, nDim, nNodes >;

    using Self = MarmotLagrangeCell< nDim, nNodes >;

    MarmotLagrangeCell() = default;
    explicit MarmotLagrangeCell( const Mat& nodes ) : nodes_( nodes ) {}
    explicit MarmotLagrangeCell( const Scalar* nodesData, int nNodesData )
    {
      if ( nNodesData != nDim * nNodes )
        throw std::invalid_argument( "MarmotLagrangeCell: invalid nodes data size" );
      nodes_ = Eigen::Map< const Mat >( nodesData );
    }

    const Mat& nodes() const { return nodes_; }
    Mat&       nodes() { return nodes_; }

    // ------------------------------------------------------------------------
    // Element-specific: shape functions and derivatives (must be specialized)
    // ------------------------------------------------------------------------
    static RowN  shapeFunctions( const Vec& s );
    static dNMat dShapeFunctions( const Vec& s );

    static RowN interpolationOperator( const Vec& s ) { return shapeFunctions( s ); }

    // ------------------------------------------------------------------------
    // Mapping and Jacobian
    // ------------------------------------------------------------------------
    Vec mapToPhysical( const Vec& s ) const { return nodes_ * shapeFunctions( s ).transpose(); }

    Eigen::Matrix< Scalar, nDim, nDim > jacobian( const Vec& s ) const
    {
      dNMat dN = dShapeFunctions( s );
      return nodes_ * dN.transpose();
    }

    // ------------------------------------------------------------------------
    // Generic exact centroid via Gauss integration
    // ------------------------------------------------------------------------
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
        Eigen::Matrix< Scalar, nDim, nDim > J    = nodes_ * dN.transpose();
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
    Scalar volume() const
    {
      Scalar V   = 0.0;
      auto   gps = GaussRule< nDim >::points();
      auto   wgs = GaussRule< nDim >::weights();

      for ( size_t i = 0; i < gps.size(); ++i ) {
        const auto& s = gps[i];
        const auto& w = wgs[i];

        dNMat                               dN   = dShapeFunctions( s );
        Eigen::Matrix< Scalar, nDim, nDim > J    = nodes_ * dN.transpose();
        const Scalar                        detJ = std::abs( J.determinant() ) * w;
        V += detJ;
      }
      return V;
    }

    Scalar area() const
    {
      static_assert( nDim == 2, "area() only makes sense in 2D." );
      return volume();
    }

    // ------------------------------------------------------------------------
// Full second moment matrix: I_kj = ∫ d_k d_j dV
// d = x - C
// ------------------------------------------------------------------------
Eigen::Matrix<Scalar, nDim, nDim> secondMoments() const
{
    Vec C = centroid();

    Eigen::Matrix<Scalar, nDim, nDim> I;
    I.setZero();

    auto gps = GaussRule<nDim>::points();
    auto wgs = GaussRule<nDim>::weights();

    for (size_t i = 0; i < gps.size(); ++i) {
        const auto& s = gps[i];
        const auto& w = wgs[i];

        dNMat dN = dShapeFunctions(s);
        Eigen::Matrix<Scalar, nDim, nDim> J = nodes_ * dN.transpose();
        Scalar detJ = std::abs(J.determinant()) * w;

        // physical coordinates
        Vec x = mapToPhysical(s);
        Vec d = x - C; // deviation from centroid

        // Accumulate outer product: d_k d_j
        I += (d * d.transpose()) * detJ;
    }

    return I;
}


    // ------------------------------------------------------------------------
    // Update vertex coordinates
    // ------------------------------------------------------------------------
    void updateVertexCoordinates( const Mat& newNodes ) { nodes_ = newNodes; }

    // ------------------------------------------------------------------------
    // Apply deformation gradient relative to centroid:
    // x_new = C + F * (x_old - C)
    // ------------------------------------------------------------------------
    void applyDeformationGradient( const Eigen::Matrix< Scalar, nDim, nDim >& F )
    {
      Vec C = centroid();
      for ( int i = 0; i < nNodes; ++i ) {
        Vec Xrel        = nodes_.col( i ) - C;
        nodes_.col( i ) = C + F * Xrel;
      }
    }

    void applyUniformDisplacement( const Vec& dX )
    {
      for ( int i = 0; i < nNodes; ++i ) {
        nodes_.col( i ) += dX;
      }
    }

    // ------------------------------------------------------------------------
    // Faces: center coords and boundary surface vectors
    //   - must be specialized for concrete topologies
    // ------------------------------------------------------------------------
    Vec getFaceCenterCoordinates( int faceId ) const
    {
      throw std::logic_error( "getFaceCenterCoordinates not implemented for this (nDim,nNodes)" );
    }

    Vec boundarySurfaceVector( int faceId ) const
    {
      throw std::logic_error( "boundarySurfaceVector not implemented for this (nDim,nNodes)" );
    }

    // ------------------------------------------------------------------------
    // Uniform subdivision: must be specialized per topology
    // ------------------------------------------------------------------------
    std::vector< Self > uniformSubdivided( int levels ) const
    {
      throw std::logic_error( "uniformSubdivided not implemented for this (nDim,nNodes)" );
    }

    // ------------------------------------------------------------------------
    // Ensight Gold cell shape name (specialized for supported element types)
    // ------------------------------------------------------------------------
    std::string getCellShape() const
    {
      throw std::logic_error( "getCellShape() not implemented for this (nDim,nNodes)" );
    }

    int getNumberOfNodes() const { return nNodes; }

    int getNumberOfDimensions() const { return nDim; }

    int getNumberOfFaces() const
    {
      throw std::logic_error( "getNumberOfFaces() not implemented for this (nDim,nNodes)" );
    }

  private:
    Mat nodes_;
  };

  // ============================================================================
  //   Specialization: Quad4  (nDim=2, nNodes=4) – Abaqus node numbering
  // ============================================================================

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

  // Faces S1..S4 in Abaqus: S1 bottom(1-2), S2 right(2-3), S3 top(3-4), S4 left(4-1)
  template <>
  inline MarmotLagrangeCell< 2, 4 >::Vec MarmotLagrangeCell< 2, 4 >::getFaceCenterCoordinates( int faceId ) const
  {
    static const int edges[4][2] = {
      { 0, 1 }, // S1
      { 1, 2 }, // S2
      { 2, 3 }, // S3
      { 3, 0 }  // S4
    };

    if ( faceId < 1 || faceId > 4 )
      throw std::out_of_range( "Quad4 faceId must be 1..4 (Abaqus)" );

    int i0 = edges[faceId - 1][0];
    int i1 = edges[faceId - 1][1];
    return 0.5 * ( nodes_.col( i0 ) + nodes_.col( i1 ) );
  }

  template <>
  inline MarmotLagrangeCell< 2, 4 >::Vec MarmotLagrangeCell< 2, 4 >::boundarySurfaceVector( int faceId ) const
  {
    static const int edges[4][2] = {
      { 0, 1 }, // S1
      { 1, 2 }, // S2
      { 2, 3 }, // S3
      { 3, 0 }  // S4
    };

    if ( faceId < 1 || faceId > 4 )
      throw std::out_of_range( "Quad4 faceId must be 1..4 (Abaqus)" );

    int i0 = edges[faceId - 1][0];
    int i1 = edges[faceId - 1][1];

    Vec e = nodes_.col( i1 ) - nodes_.col( i0 );
    return Vec( e[1], -e[0] ); // outward normal * edge length (for CCW nodes)
  }

  template <>
  inline std::vector< MarmotLagrangeCell< 2, 4 > > MarmotLagrangeCell< 2, 4 >::uniformSubdivided( int levels ) const
  {
    std::vector< Self > elems = { *this };

    for ( int l = 0; l < levels; ++l ) {
      std::vector< Self > next;
      next.reserve( elems.size() * 4 );

      for ( const auto& e : elems ) {
        Vec s;
        Vec grid[3][3];

        for ( int j = 0; j < 3; ++j )
          for ( int i = 0; i < 3; ++i ) {
            s << -1.0 + i, -1.0 + j;
            grid[i][j] = e.mapToPhysical( s );
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
      }

      elems.swap( next );
    }

    return elems;
  }

  // ============================================================================
  //   Specialization: Hex8  (nDim=3, nNodes=8) – Abaqus node & face numbering
  // ============================================================================

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

  // Abaqus face numbering:
  // S1: 1,2,3,4  bottom
  // S2: 5,6,7,8  top
  // S3: 1,2,6,5  front
  // S4: 2,3,7,6  right
  // S5: 3,4,8,7  back
  // S6: 4,1,5,8  left
  template <>
  inline MarmotLagrangeCell< 3, 8 >::Vec MarmotLagrangeCell< 3, 8 >::getFaceCenterCoordinates( int faceId ) const
  {
    static const int faces[6][4] = {
      { 0, 1, 2, 3 }, // S1
      { 4, 5, 6, 7 }, // S2
      { 0, 1, 5, 4 }, // S3
      { 1, 2, 6, 5 }, // S4
      { 2, 3, 7, 6 }, // S5
      { 3, 0, 4, 7 }  // S6
    };

    if ( faceId < 1 || faceId > 6 )
      throw std::out_of_range( "Hex8 faceId must be 1..6 (Abaqus)" );

    Vec c = Vec::Zero();
    for ( int i = 0; i < 4; ++i )
      c += nodes_.col( faces[faceId - 1][i] );

    return c / 4.0;
  }

  template <>
  inline MarmotLagrangeCell< 3, 8 >::Vec MarmotLagrangeCell< 3, 8 >::boundarySurfaceVector( int faceId ) const
  {
    static const int faces[6][4] = {
      { 0, 1, 2, 3 }, // S1
      { 4, 5, 6, 7 }, // S2
      { 0, 1, 5, 4 }, // S3
      { 1, 2, 6, 5 }, // S4
      { 2, 3, 7, 6 }, // S5
      { 3, 0, 4, 7 }  // S6
    };

    if ( faceId < 1 || faceId > 6 )
      throw std::out_of_range( "Hex8 faceId must be 1..6 (Abaqus)" );

    int a = faces[faceId - 1][0];
    int b = faces[faceId - 1][1];
    int c = faces[faceId - 1][2];
    int d = faces[faceId - 1][3];

    Vec p1 = nodes_.col( a );
    Vec p2 = nodes_.col( b );
    Vec p3 = nodes_.col( c );
    Vec p4 = nodes_.col( d );

    Vec n1 = ( p2 - p1 ).cross( p3 - p1 );
    Vec n2 = ( p4 - p1 ).cross( p3 - p1 );

    return 0.5 * ( n1 + n2 );
  }

  template <>
  inline std::vector< MarmotLagrangeCell< 3, 8 > > MarmotLagrangeCell< 3, 8 >::uniformSubdivided( int levels ) const
  {
    std::vector< Self > elems = { *this };

    for ( int l = 0; l < levels; ++l ) {
      std::vector< Self > next;
      next.reserve( elems.size() * 8 );

      for ( const auto& e : elems ) {
        Vec s;
        Vec grid[3][3][3];

        for ( int k = 0; k < 3; ++k )
          for ( int j = 0; j < 3; ++j )
            for ( int i = 0; i < 3; ++i ) {
              s << -1.0 + i, -1.0 + j, -1.0 + k;
              grid[i][j][k] = e.mapToPhysical( s );
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
      }

      elems.swap( next );
    }

    return elems;
  }

  // ============================================================================
  // Ensight Gold cell shape specialization: Quad4
  // ============================================================================
  template <>
  inline std::string MarmotLagrangeCell< 2, 4 >::getCellShape() const
  {
    return "quad4"; // Ensight Gold notation
  }

  // ============================================================================
  // Ensight Gold cell shape specialization: Hex8
  // ============================================================================
  template <>
  inline std::string MarmotLagrangeCell< 3, 8 >::getCellShape() const
  {
    return "hexa8"; // Ensight Gold notation
  }

  // number of faces:
  template <>
  inline int MarmotLagrangeCell< 2, 4 >::getNumberOfFaces() const
  {
    return 4;
  }

  template <>
  inline int MarmotLagrangeCell< 3, 8 >::getNumberOfFaces() const
  {
    return 6;
  }

  // ============================================================================
  // Convenient aliases if you still want Quad4 / Hex8 names
  // ============================================================================

  using Quad4 = MarmotLagrangeCell< 2, 4 >;
  using Hex8  = MarmotLagrangeCell< 3, 8 >;

} // namespace Marmot::Meshfree
