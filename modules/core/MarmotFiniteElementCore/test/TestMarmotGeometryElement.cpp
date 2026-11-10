#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::Math::Differentiation::Numerical;

// ---------------------------------------------------------------------------------------------
// Tetra4 / Tetra10: neither is used by any element tested elsewhere in this codebase, so their
// shape functions have never been exercised. Rather than trust a hand-transcribed formula (this
// file's Tetra10 dNdXi in particular is a long hand-expanded product-rule derivative, exactly the
// kind of code most prone to a transcription slip), correctness is checked two ways that need no
// prior knowledge of nodal coordinates: N must partition unity (and reduce to the Kronecker delta
// at each of Tetra4's easily-derived nodal points), and dNdXi must equal the numerical derivative
// of N with respect to the same natural coordinates.
// ---------------------------------------------------------------------------------------------

void testTetra4PartitionOfUnityAndKroneckerDeltaAtNodes()
{
  using namespace Marmot::Discretizations::Elements::Spatial3D;

  // Derived directly from N: N0=1-xi0-xi1-xi2, N1=xi0, N2=xi1, N3=xi2.
  const std::vector< Eigen::Vector3d > nodalXi = { { 0.0, 0.0, 0.0 },
                                                   { 1.0, 0.0, 0.0 },
                                                   { 0.0, 1.0, 0.0 },
                                                   { 0.0, 0.0, 1.0 } };

  for ( int i = 0; i < 4; i++ ) {
    const auto N = Tetra4::N( nodalXi[i] );
    throwExceptionOnFailure( checkIfEqual( N.sum(), 1.0, 1e-12 ), "Tetra4::N() does not partition unity." );
    for ( int j = 0; j < 4; j++ ) {
      const double expected = ( i == j ) ? 1.0 : 0.0;
      throwExceptionOnFailure( checkIfEqual( N( j ), expected, 1e-12 ),
                               "Tetra4::N() does not reduce to the Kronecker delta at its own nodes." );
    }
  }
}

void testTetra4DNdXiMatchesNumericalDifferentiationOfN()
{
  using namespace Marmot::Discretizations::Elements::Spatial3D;

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );

  NumDiff::vector_to_vector_function_type Nfunc = []( const Eigen::VectorXd& xiVec ) -> Eigen::VectorXd {
    return Tetra4::N( Eigen::Vector3d( xiVec ) );
  };

  const Eigen::MatrixXd numJacobian = NumDiff::centralDifference( Nfunc, xi ); // (dN_i/dXi_j), 4x3

  const auto analyticDNdXi = Tetra4::dNdXi( xi ); // (dN_i/dXi_j) stored as (nDim x nNodes) = dNdXi(j,i)

  throwExceptionOnFailure( checkIfEqual( numJacobian, Eigen::MatrixXd( analyticDNdXi.transpose() ), 1e-6 ),
                           "Tetra4::dNdXi() does not match the numerical derivative of Tetra4::N()." );
}

void testTetra10PartitionOfUnityAtSamplePoints()
{
  using namespace Marmot::Discretizations::Elements::Spatial3D;

  const std::vector< Eigen::Vector3d > samples = { { 0.0, 0.0, 0.0 },
                                                   { 1.0, 0.0, 0.0 },
                                                   { 0.0, 1.0, 0.0 },
                                                   { 0.0, 0.0, 1.0 },
                                                   { 0.2, 0.3, 0.1 },
                                                   { 0.25, 0.25, 0.25 } };

  for ( const auto& xi : samples ) {
    const auto N = Tetra10::N( xi );
    throwExceptionOnFailure( checkIfEqual( N.sum(), 1.0, 1e-10 ), "Tetra10::N() does not partition unity." );
  }
}

void testTetra10DNdXiMatchesNumericalDifferentiationOfN()
{
  using namespace Marmot::Discretizations::Elements::Spatial3D;

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );

  NumDiff::vector_to_vector_function_type Nfunc = []( const Eigen::VectorXd& xiVec ) -> Eigen::VectorXd {
    return Tetra10::N( Eigen::Vector3d( xiVec ) );
  };

  const Eigen::MatrixXd numJacobian   = NumDiff::centralDifference( Nfunc, xi );
  const auto            analyticDNdXi = Tetra10::dNdXi( xi );

  throwExceptionOnFailure( checkIfEqual( numJacobian, Eigen::MatrixXd( analyticDNdXi.transpose() ), 1e-6 ),
                           "Tetra10::dNdXi() does not match the numerical derivative of Tetra10::N() -- likely "
                           "a transcription error in the hand-expanded derivative." );
}

// ---------------------------------------------------------------------------------------------
// BGreen(dNdX, F): at F = Identity, the Green-Lagrange strain operator reduces algebraically to
// the standard linear B-operator (every cross term multiplies a zero off-diagonal entry of F, and
// every retained term multiplies F's unit diagonal) -- an exact identity, not an approximation,
// and one that exercises Marmot::Discretizations::Elements::GeometryElement's BGreen() specialization (which has no
// callers anywhere in the codebase) for every shape.
// ---------------------------------------------------------------------------------------------

template < int nDim, int nNodes >
void checkBGreenAtIdentityMatchesB( const std::vector< double >&            nodeCoordsVec,
                                    const Eigen::Matrix< double, nDim, 1 >& xi,
                                    const std::string&                      label )
{
  Marmot::Discretizations::Elements::GeometryElement< nDim, nNodes > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const auto dNdXi = geo.dNdXi( xi );
  const auto J     = geo.Jacobian( dNdXi );
  const auto dNdX  = geo.dNdX( dNdXi, J.inverse() );

  const auto B      = geo.B( dNdX );
  const auto BGreen = geo.BGreen( dNdX, Eigen::Matrix< double, nDim, nDim >::Identity() );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( BGreen ), Eigen::MatrixXd( B ), 1e-10 ),
                           label + ": BGreen(dNdX, Identity) must equal B(dNdX)." );
}

void testBGreenAtIdentityMatchesBForEveryShape()
{
  const std::vector< double > quad4Coords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 2, 4 >( quad4Coords, Eigen::Vector2d( 0.2, 0.3 ), "Quad4" );

  const std::vector< double > quad8Coords =
    { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5 };
  checkBGreenAtIdentityMatchesB< 2, 8 >( quad8Coords, Eigen::Vector2d( 0.2, 0.3 ), "Quad8" );

  const std::vector< double > hexa8Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 3, 8 >( hexa8Coords, Eigen::Vector3d( 0.2, 0.3, 0.1 ), "Hexa8" );

  const std::vector< double > hexa20Coords = {
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, // bottom corners
    0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, // top corners
    0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5, 0.0, // bottom-face edge midsides
    0.5, 0.0, 1.0, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0, 0.0, 0.5, 1.0, // top-face edge midsides
    0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5  // vertical edge midsides
  };
  checkBGreenAtIdentityMatchesB< 3, 20 >( hexa20Coords, Eigen::Vector3d( 0.2, 0.3, 0.1 ), "Hexa20" );

  const std::vector< double > tetra4Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  checkBGreenAtIdentityMatchesB< 3, 4 >( tetra4Coords, Eigen::Vector3d( 0.2, 0.2, 0.2 ), "Tetra4" );
}

// ---------------------------------------------------------------------------------------------
// B_bar(dNdX, dNdX0): when evaluated at the SAME point (dNdX0 == dNdX), the B-bar correction terms
// (B4, B6, B8) vanish identically and B_bar algebraically collapses to the standard B operator --
// again an exact identity, not an approximation.
// ---------------------------------------------------------------------------------------------

void testBBarAtEqualPointsMatchesBForHexa8()
{
  const std::vector< double > hexa8Coords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  Marmot::Discretizations::Elements::GeometryElement< 3, 8 > geo;
  geo.assignNodeCoordinates( hexa8Coords.data() );

  const Eigen::Vector3d xi( 0.2, 0.3, 0.1 );
  const auto            dNdXi = geo.dNdXi( xi );
  const auto            J     = geo.Jacobian( dNdXi );
  const auto            dNdX  = geo.dNdX( dNdXi, J.inverse() );

  const auto B     = geo.B( dNdX );
  const auto B_bar = geo.B_bar( dNdX, dNdX );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( B_bar ), Eigen::MatrixXd( B ), 1e-10 ),
                           "B_bar(dNdX, dNdX) must equal B(dNdX) when evaluated at the same point." );
}

// ---------------------------------------------------------------------------------------------
// B_axisymmetric: no algebraic identity collapses this to the standard B operator (it introduces
// a genuinely new N/r hoop-strain term), so it is checked against hand-computed values instead.
// ---------------------------------------------------------------------------------------------

void testBAxisymmetricMatchesHandDerivedValuesQuad4()
{
  Marmot::Discretizations::Elements::GeometryElement< 2, 4 >::dNdXiSized dNdX;
  dNdX << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0; // arbitrary but fixed 2x4 values

  Marmot::Discretizations::Elements::GeometryElement< 2, 4 >::NSized N;
  N << 0.1, 0.2, 0.3, 0.4;

  const Eigen::Vector2d x_gauss( 2.0, 0.0 ); // r = 2

  Marmot::Discretizations::Elements::GeometryElement< 2, 4 > geo;
  const auto                                                 B = geo.B_axisymmetric( dNdX, N, x_gauss );

  throwExceptionOnFailure( B.rows() == 4 && B.cols() == 8, "B_axisymmetric() returned the wrong size." );

  for ( int i = 0; i < 4; i++ ) {
    throwExceptionOnFailure( checkIfEqual( B( 0, 2 * i ), dNdX( 0, i ), 1e-12 ), "B_axisymmetric row 0 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 0, 2 * i + 1 ), 0.0, 1e-12 ), "B_axisymmetric row 0 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 1, 2 * i + 1 ), dNdX( 1, i ), 1e-12 ), "B_axisymmetric row 1 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 1, 2 * i ), 0.0, 1e-12 ), "B_axisymmetric row 1 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 2, 2 * i ), N( i ) / x_gauss( 0 ), 1e-12 ),
                             "B_axisymmetric row 2 (hoop strain N/r) mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 3, 2 * i ), dNdX( 1, i ), 1e-12 ), "B_axisymmetric row 3 mismatch." );
    throwExceptionOnFailure( checkIfEqual( B( 3, 2 * i + 1 ), dNdX( 0, i ), 1e-12 ), "B_axisymmetric row 3 mismatch." );
  }
}

// ---------------------------------------------------------------------------------------------
// Element geometry depending on the reference domain: center, characteristic element length and
// volume-equivalent length. The cube [-1, 1]^nDim and the unit simplex differ in all three.
// ---------------------------------------------------------------------------------------------

namespace {
  using Marmot::Discretizations::Elements::GeometryElement;
  using Marmot::Discretizations::Elements::ReferenceDomain;

  static_assert( GeometryElement< 2, 3 >::referenceDomain == ReferenceDomain::Triangle, "Tri3 is a triangle." );
  static_assert( GeometryElement< 3, 4 >::referenceDomain == ReferenceDomain::Tetrahedron, "Tetra4 is a tetrahedron." );
  static_assert( GeometryElement< 2, 8 >::referenceDomain == ReferenceDomain::Quadrilateral, "Quad8 is a quad." );
  static_assert( GeometryElement< 3, 20 >::referenceDomain == ReferenceDomain::Hexahedron, "Hexa20 is a hexahedron." );
  static_assert( GeometryElement< 3, 10 >::interpolationOrder == 2, "Tetra10 is quadratic." );

  using namespace Marmot::Discretizations::Elements;
  static_assert( isSimplex( ReferenceDomain::Triangle ) && isSimplex( ReferenceDomain::Tetrahedron ) &&
                   !isSimplex( ReferenceDomain::Line ) && !isSimplex( ReferenceDomain::Quadrilateral ) &&
                   !isSimplex( ReferenceDomain::Hexahedron ),
                 "Only the triangle and the tetrahedron are simplices." );
  static_assert( getReferenceMeasure( ReferenceDomain::Line ) == 2.0 &&
                   getReferenceMeasure( ReferenceDomain::Triangle ) == 1.0 / 2.0 &&
                   getReferenceMeasure( ReferenceDomain::Quadrilateral ) == 4.0 &&
                   getReferenceMeasure( ReferenceDomain::Tetrahedron ) == 1.0 / 6.0 &&
                   getReferenceMeasure( ReferenceDomain::Hexahedron ) == 8.0,
                 "Reference measures: 2^nDim for the cube, 1 / nDim! for the unit simplex." );
  static_assert( getReferenceCentroidCoordinate( ReferenceDomain::Quadrilateral ) == 0.0 &&
                   getReferenceCentroidCoordinate( ReferenceDomain::Hexahedron ) == 0.0 &&
                   getReferenceCentroidCoordinate( ReferenceDomain::Triangle ) == 1.0 / 3.0 &&
                   getReferenceCentroidCoordinate( ReferenceDomain::Tetrahedron ) == 1.0 / 4.0,
                 "Reference centroids: the origin of the cube, 1 / ( nDim + 1 ) for the unit simplex." );

  /// distance of point @p p from the line (2D) or plane (3D) through the points @p face
  template < int nDim >
  double distanceFromFace( const Eigen::Matrix< double, nDim, 1 >&                p,
                           const std::vector< Eigen::Matrix< double, nDim, 1 > >& face )
  {
    if constexpr ( nDim == 2 ) {
      const Eigen::Vector2d t = face[1] - face[0], d = p - face[0];
      return std::abs( t( 0 ) * d( 1 ) - t( 1 ) * d( 0 ) ) / t.norm();
    }
    else {
      const Eigen::Vector3d n = ( face[1] - face[0] ).cross( face[2] - face[0] );
      return std::abs( n.dot( p - face[0] ) ) / n.norm();
    }
  }
} // namespace

template < int nDim >
void checkGeometryOfLinearSimplex( const std::vector< double >& coordinates, double measure )
{
  constexpr int nNodes = nDim + 1;
  using Simplex        = GeometryElement< nDim, nNodes >;
  using Point          = Eigen::Matrix< double, nDim, 1 >;

  double firstNumbering = 0.0;
  for ( int shift = 0; shift < nNodes; shift++ ) {
    // a cyclic shift of the four nodes of a tetrahedron is an odd permutation, which would invert it; swapping two
    // nodes as well keeps the orientation
    std::vector< int > permutation( nNodes );
    for ( int i = 0; i < nNodes; i++ )
      permutation[i] = ( i + shift ) % nNodes;
    if constexpr ( nNodes % 2 == 0 )
      if ( shift % 2 == 1 )
        std::swap( permutation[1], permutation[2] );

    std::vector< double > shifted( coordinates.size() );
    for ( int i = 0; i < nNodes; i++ )
      for ( int d = 0; d < nDim; d++ )
        shifted[i * nDim + d] = coordinates[permutation[i] * nDim + d];

    Simplex simplex;
    simplex.assignNodeCoordinates( shifted.data() );

    std::vector< Point > nodes( nNodes );
    for ( int i = 0; i < nNodes; i++ )
      nodes[i] = Eigen::Map< const Point >( &shifted[i * nDim] );

    // the center is the centroid, the mean of the nodes, not node 0 at the origin of the unit simplex
    Point centroid = Point::Zero();
    for ( const auto& node : nodes )
      centroid += node / nNodes;
    throwExceptionOnFailure( ( simplex.center() - centroid ).norm() < 1e-14,
                             "The center of a linear simplex must be its centroid." );

    // the characteristic element length is 1 / sqrt( sum_i h_i^-2 ), with the altitude h_i opposite node i
    double sumOfInverseSquaredAltitudes = 0.0;
    for ( int i = 0; i < nNodes; i++ ) {
      std::vector< Point > face;
      for ( int j = 0; j < nNodes; j++ )
        if ( j != i )
          face.push_back( nodes[j] );
      sumOfInverseSquaredAltitudes += std::pow( distanceFromFace< nDim >( nodes[i], face ), -2 );
    }
    const double length = simplex.characteristicLengthForCriticalTimeStepAt( Simplex::referenceCentroid() );
    throwExceptionOnFailure( checkIfEqual( length, 1.0 / std::sqrt( sumOfInverseSquaredAltitudes ), 1e-12 ),
                             "The characteristic element length of a linear simplex must be 1 / sqrt( sum h^-2 )." );
    if ( shift == 0 )
      firstNumbering = length;
    throwExceptionOnFailure( std::abs( length - firstNumbering ) < 1e-12 * firstNumbering,
                             "The characteristic element length must not depend on the node numbering." );

    // the volume-equivalent length is sqrt( A ) or cbrt( V )
    const double detJ = simplex.Jacobian( simplex.dNdXi( Simplex::referenceCentroid() ) ).determinant();
    throwExceptionOnFailure( checkIfEqual( Simplex::volumeEquivalentLength( detJ ),
                                           nDim == 2 ? std::sqrt( measure ) : std::cbrt( measure ),
                                           1e-12 ),
                             "The volume-equivalent length of a linear simplex must be sqrt( A ) or cbrt( V )." );
  }
}

void testGeometryOfLinearSimplices()
{
  checkGeometryOfLinearSimplex< 2 >( { 0, 0, 2, 0, 0.5, 1 }, 1.0 );
  checkGeometryOfLinearSimplex< 2 >( { 0, 0, 1, 0, 0.5, std::sqrt( 3. ) / 2 }, std::sqrt( 3. ) / 4 );
  checkGeometryOfLinearSimplex< 3 >( { 0, 0, 0, 2, 0, 0, 0.5, 1, 0, 0, 0.5, 3 }, 1.0 );
  checkGeometryOfLinearSimplex< 3 >( { 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1 }, 8.0 / 3.0 );
}

void testGeometryOfQuadrilateralsAndHexahedra()
{
  // a 4 x 0.5 rectangle: the characteristic element length is the smallest extent, 0.5
  GeometryElement< 2, 4 >     quad;
  const std::vector< double > quadCoordinates = { 0, 0, 4, 0, 4, 0.5, 0, 0.5 };
  quad.assignNodeCoordinates( quadCoordinates.data() );
  const Eigen::Vector2d origin = Eigen::Vector2d::Zero();
  throwExceptionOnFailure( ( quad.center() - Eigen::Vector2d( 2, 0.25 ) ).norm() < 1e-14,
                           "The center of a quadrilateral must be the image of the origin." );
  throwExceptionOnFailure( checkIfEqual( quad.characteristicLengthForCriticalTimeStepAt( origin ), 0.5, 1e-14 ),
                           "The characteristic element length of a rectangle must be its smaller side." );
  const double quadDetJ = quad.Jacobian( quad.dNdXi( origin ) ).determinant();
  throwExceptionOnFailure( checkIfEqual( GeometryElement< 2, 4 >::volumeEquivalentLength( quadDetJ ),
                                         std::sqrt( 2.0 ),
                                         1e-14 ),
                           "The volume-equivalent length of a quadrilateral must be sqrt( A )." );

  // the unit cube: all three lengths are 1
  GeometryElement< 3, 8 >     hexa;
  const std::vector< double > hexaCoordinates = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                                  0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
  hexa.assignNodeCoordinates( hexaCoordinates.data() );
  const Eigen::Vector3d origin3D = Eigen::Vector3d::Zero();
  throwExceptionOnFailure( ( hexa.center() - Eigen::Vector3d::Constant( 0.5 ) ).norm() < 1e-14,
                           "The center of the unit cube must be ( 0.5, 0.5, 0.5 )." );
  throwExceptionOnFailure( checkIfEqual( hexa.characteristicLengthForCriticalTimeStepAt( origin3D ), 1.0, 1e-14 ),
                           "The characteristic element length of the unit cube must be 1." );
  const double hexaDetJ = hexa.Jacobian( hexa.dNdXi( origin3D ) ).determinant();
  throwExceptionOnFailure( checkIfEqual( GeometryElement< 3, 8 >::volumeEquivalentLength( hexaDetJ ), 1.0, 1e-14 ),
                           "The volume-equivalent length of the unit cube must be 1." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testTetra4PartitionOfUnityAndKroneckerDeltaAtNodes,
    testTetra4DNdXiMatchesNumericalDifferentiationOfN,
    testTetra10PartitionOfUnityAtSamplePoints,
    testTetra10DNdXiMatchesNumericalDifferentiationOfN,
    testBGreenAtIdentityMatchesBForEveryShape,
    testBBarAtEqualPointsMatchesBForHexa8,
    testBAxisymmetricMatchesHandDerivedValuesQuad4,
    testGeometryOfLinearSimplices,
    testGeometryOfQuadrilateralsAndHexahedra,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
