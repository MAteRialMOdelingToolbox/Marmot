#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement;

// ---------------------------------------------------------------------------------------------
// BoundaryElement has no public accessors for its (private) condensed boundary coordinates or its
// parent<->boundary index maps, so the tests below never rely on them directly. Instead, expected
// values are derived independently from the KNOWN parent element geometry (a unit square / unit
// cube) and cross-checked against what assembleIntoParent*() actually writes into a parent-sized
// vector -- which is the same public interface every element in this codebase uses.
// ---------------------------------------------------------------------------------------------

namespace {

  const std::vector< double > unitSquareCoords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  const std::vector< double > unitCubeCoords = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  // Assembles a boundary-sized vectorial quantity into a zeroed parent-sized vector and reports
  // which parent nodes ended up with a nonzero entry (i.e., which nodes lie on this face).
  struct AssembledVectorial {
    Eigen::VectorXd    parentVector;
    std::vector< int > activeNodes;
  };

  AssembledVectorial assembleVectorial( const Eigen::VectorXd& boundaryVector,
                                        BoundaryElement&       boundaryEl,
                                        int                    nNodesParent,
                                        int                    nDim )
  {
    AssembledVectorial result;
    result.parentVector = Eigen::VectorXd::Zero( nNodesParent * nDim );
    boundaryEl.assembleIntoParentVectorial( boundaryVector, result.parentVector );

    for ( int n = 0; n < nNodesParent; n++ ) {
      bool nonzero = false;
      for ( int d = 0; d < nDim; d++ )
        if ( std::abs( result.parentVector( n * nDim + d ) ) > 1e-10 )
          nonzero = true;
      if ( nonzero )
        result.activeNodes.push_back( n );
    }
    return result;
  }

} // namespace

// ---------------------------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------------------------

void testConstructorThrowsForUnsupportedParentShape()
{
  const Eigen::VectorXd coords = Eigen::VectorXd::Zero( 12 ); // Tetra4: 4 nodes * 3 dim

  bool threw = false;
  try {
    BoundaryElement boundaryEl( Tetra4, 0, 3, coords );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "BoundaryElement must throw for an unsupported parent shape." );
}

// ---------------------------------------------------------------------------------------------
// computeScalarLoadVector(): integrates a unit scalar field, so its total (summed via
// assembleIntoParentScalar) must equal the physical length/area of the face, by partition of
// unity (sum_i N_i = 1 everywhere).
// ---------------------------------------------------------------------------------------------

void testComputeScalarLoadVectorIntegratesToFaceLengthQuad4()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );

  for ( int face = 1; face <= 4; face++ ) {
    BoundaryElement       boundaryEl( Quad4, face, 2, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeScalarLoadVector();

    Eigen::VectorXd parentScalar = Eigen::VectorXd::Zero( 4 ); // 4 parent nodes
    boundaryEl.assembleIntoParentScalar( Pk, parentScalar );

    throwExceptionOnFailure( checkIfEqual( parentScalar.sum(), 1.0, 1e-10 ),
                             "computeScalarLoadVector() does not integrate to the unit square's edge "
                             "length (face " +
                               std::to_string( face ) + ")." );
  }
}

void testComputeScalarLoadVectorIntegratesToFaceAreaHexa8()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitCubeCoords.data(), 24 );

  for ( int face = 1; face <= 6; face++ ) {
    BoundaryElement       boundaryEl( Hexa8, face, 3, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeScalarLoadVector();

    Eigen::VectorXd parentScalar = Eigen::VectorXd::Zero( 8 ); // 8 parent nodes
    boundaryEl.assembleIntoParentScalar( Pk, parentScalar );

    throwExceptionOnFailure( checkIfEqual( parentScalar.sum(), 1.0, 1e-10 ),
                             "computeScalarLoadVector() does not integrate to the unit cube's face "
                             "area (face " +
                               std::to_string( face ) + ")." );
  }
}

// ---------------------------------------------------------------------------------------------
// computeVectorialLoadVector(): integrates a constant traction vector, so the resultant force
// (summed per direction) must equal direction * face length/area.
// ---------------------------------------------------------------------------------------------

void testComputeVectorialLoadVectorIntegratesToDirectionTimesLengthQuad4()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  const Eigen::Vector2d                     direction( 2.0, -3.0 );

  for ( int face = 1; face <= 4; face++ ) {
    BoundaryElement       boundaryEl( Quad4, face, 2, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeVectorialLoadVector( direction );

    const auto assembled = assembleVectorial( Pk, boundaryEl, 4, 2 );
    throwExceptionOnFailure( assembled.activeNodes.size() == 2,
                             "Exactly 2 parent nodes should lie on a Quad4 edge (face " + std::to_string( face ) +
                               ")." );

    Eigen::Vector2d resultant = Eigen::Vector2d::Zero();
    for ( int n : assembled.activeNodes )
      resultant += assembled.parentVector.segment< 2 >( n * 2 );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( resultant ), Eigen::MatrixXd( direction ), 1e-10 ),
                             "computeVectorialLoadVector() resultant does not match direction * edge length "
                             "(face " +
                               std::to_string( face ) + ")." );
  }
}

void testComputeVectorialLoadVectorIntegratesToDirectionTimesAreaHexa8()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitCubeCoords.data(), 24 );
  const Eigen::Vector3d                     direction( 1.0, 2.0, -1.5 );

  for ( int face = 1; face <= 6; face++ ) {
    BoundaryElement       boundaryEl( Hexa8, face, 3, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeVectorialLoadVector( direction );

    const auto assembled = assembleVectorial( Pk, boundaryEl, 8, 3 );
    throwExceptionOnFailure( assembled.activeNodes.size() == 4,
                             "Exactly 4 parent nodes should lie on a Hexa8 face (face " + std::to_string( face ) +
                               ")." );

    Eigen::Vector3d resultant = Eigen::Vector3d::Zero();
    for ( int n : assembled.activeNodes )
      resultant += assembled.parentVector.segment< 3 >( n * 3 );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( resultant ), Eigen::MatrixXd( direction ), 1e-10 ),
                             "computeVectorialLoadVector() resultant does not match direction * face area "
                             "(face " +
                               std::to_string( face ) + ")." );
  }
}

// ---------------------------------------------------------------------------------------------
// computeSurfaceNormalVectorialLoadVector(): the resultant must point along the OUTWARD normal
// (determined independently here from which way the face centroid lies relative to the parent
// element's centroid) with magnitude equal to the face length/area.
// ---------------------------------------------------------------------------------------------

void testComputeSurfaceNormalVectorialLoadVectorPointsOutwardQuad4()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  const Eigen::Vector2d                     parentCentroid( 0.5, 0.5 );

  for ( int face = 1; face <= 4; face++ ) {
    BoundaryElement       boundaryEl( Quad4, face, 2, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeSurfaceNormalVectorialLoadVector();

    const auto assembled = assembleVectorial( Pk, boundaryEl, 4, 2 );
    throwExceptionOnFailure( assembled.activeNodes.size() == 2,
                             "Exactly 2 parent nodes should lie on a Quad4 edge (face " + std::to_string( face ) +
                               ")." );

    const Eigen::Vector2d P0( unitSquareCoords[assembled.activeNodes[0] * 2],
                              unitSquareCoords[assembled.activeNodes[0] * 2 + 1] );
    const Eigen::Vector2d P1( unitSquareCoords[assembled.activeNodes[1] * 2],
                              unitSquareCoords[assembled.activeNodes[1] * 2 + 1] );

    const Eigen::Vector2d tangent      = P1 - P0;
    const double          edgeLength   = tangent.norm();
    Eigen::Vector2d       candidate    = Eigen::Vector2d( tangent( 1 ), -tangent( 0 ) ).normalized();
    const Eigen::Vector2d edgeMidpoint = 0.5 * ( P0 + P1 );
    if ( ( edgeMidpoint - parentCentroid ).dot( candidate ) < 0.0 )
      candidate = -candidate;

    const Eigen::Vector2d expected = candidate * edgeLength;

    Eigen::Vector2d resultant = Eigen::Vector2d::Zero();
    for ( int n : assembled.activeNodes )
      resultant += assembled.parentVector.segment< 2 >( n * 2 );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( resultant ), Eigen::MatrixXd( expected ), 1e-10 ),
                             "computeSurfaceNormalVectorialLoadVector() does not point along the outward "
                             "normal with the correct magnitude (face " +
                               std::to_string( face ) + ")." );
  }
}

void testComputeSurfaceNormalVectorialLoadVectorPointsOutwardHexa8()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitCubeCoords.data(), 24 );
  const Eigen::Vector3d                     parentCentroid( 0.5, 0.5, 0.5 );

  for ( int face = 1; face <= 6; face++ ) {
    BoundaryElement       boundaryEl( Hexa8, face, 3, parentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeSurfaceNormalVectorialLoadVector();

    const auto assembled = assembleVectorial( Pk, boundaryEl, 8, 3 );
    throwExceptionOnFailure( assembled.activeNodes.size() == 4,
                             "Exactly 4 parent nodes should lie on a Hexa8 face (face " + std::to_string( face ) +
                               ")." );

    std::vector< Eigen::Vector3d > corners;
    for ( int n : assembled.activeNodes )
      corners.emplace_back( unitCubeCoords[n * 3], unitCubeCoords[n * 3 + 1], unitCubeCoords[n * 3 + 2] );

    Eigen::Vector3d       candidate    = ( corners[1] - corners[0] ).cross( corners[2] - corners[0] ).normalized();
    const Eigen::Vector3d faceCentroid = 0.25 * ( corners[0] + corners[1] + corners[2] + corners[3] );
    if ( ( faceCentroid - parentCentroid ).dot( candidate ) < 0.0 )
      candidate = -candidate;

    const double          unitCubeFaceArea = 1.0;
    const Eigen::Vector3d expected         = candidate * unitCubeFaceArea;

    Eigen::Vector3d resultant = Eigen::Vector3d::Zero();
    for ( int n : assembled.activeNodes )
      resultant += assembled.parentVector.segment< 3 >( n * 3 );

    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( resultant ), Eigen::MatrixXd( expected ), 1e-10 ),
                             "computeSurfaceNormalVectorialLoadVector() does not point along the outward "
                             "normal with the correct magnitude (face " +
                               std::to_string( face ) + ")." );
  }
}

// ---------------------------------------------------------------------------------------------
// computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements(): for a face at CONSTANT radial
// coordinate r (so the 2*pi*r factor is the same at every quadrature point), the result must be
// exactly the plain surface-normal load scaled by 2*pi*r.
// ---------------------------------------------------------------------------------------------

void testComputeSurfaceNormalVectorialLoadVectorForAxisymmetricElementsScalesByTwoPiR()
{
  // A square with r in [1, 2], z in [0, 1]: face 1 (r=1 -> r=2 at z=0, from the node ordering
  // below) is not constant-r, but the LEFT edge (nodes 3-0, r=1) and RIGHT edge (nodes 1-2, r=2)
  // are. We scan all faces and only assert on the ones with constant r.
  const std::vector< double >               coordsVec = { 1.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 1.0 };
  const Eigen::Map< const Eigen::VectorXd > parentCoords( coordsVec.data(), 8 );

  bool checkedAtLeastOneFace = false;

  for ( int face = 1; face <= 4; face++ ) {
    BoundaryElement       boundaryEl( Quad4, face, 2, parentCoords );
    const Eigen::VectorXd plain        = boundaryEl.computeSurfaceNormalVectorialLoadVector();
    const Eigen::VectorXd axisymmetric = boundaryEl.computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements();

    const auto   assembled = assembleVectorial( plain, boundaryEl, 4, 2 );
    const double r0        = coordsVec[assembled.activeNodes[0] * 2];
    const double r1        = coordsVec[assembled.activeNodes[1] * 2];

    if ( std::abs( r0 - r1 ) > 1e-12 )
      continue; // r varies along this face; the simple scaling relation does not hold here

    checkedAtLeastOneFace = true;
    const double factor   = 2.0 * Constants::Pi * r0;
    throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( axisymmetric ), Eigen::MatrixXd( factor * plain ), 1e-8 ),
                             "computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements() does not equal "
                             "2*pi*r times the plain surface-normal load for a constant-r face." );
  }

  throwExceptionOnFailure( checkedAtLeastOneFace, "Test setup error: no constant-r face was found to check." );
}

// ---------------------------------------------------------------------------------------------
// Quad8 / Hexa20: exercise the remaining constructor branches. Only the partition-of-unity total
// is checked here (no corner-geometry assumptions), since these boundary shapes carry midside
// nodes.
// ---------------------------------------------------------------------------------------------

void testQuad8AndHexa20BoundariesIntegrateToFaceLengthAndArea()
{
  // Quad8: unit square with midside nodes at exact edge midpoints.
  const std::vector< double > quad8Coords =
    { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5 };
  const Eigen::Map< const Eigen::VectorXd > quad8ParentCoords( quad8Coords.data(), 16 );

  for ( int face = 1; face <= 4; face++ ) {
    BoundaryElement       boundaryEl( Quad8, face, 2, quad8ParentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeScalarLoadVector();

    Eigen::VectorXd parentScalar = Eigen::VectorXd::Zero( 8 );
    boundaryEl.assembleIntoParentScalar( Pk, parentScalar );

    throwExceptionOnFailure( checkIfEqual( parentScalar.sum(), 1.0, 1e-10 ),
                             "Quad8 boundary computeScalarLoadVector() does not integrate to the edge length "
                             "(face " +
                               std::to_string( face ) + ")." );
  }

  // Hexa20: unit cube with edge-midside nodes at exact midpoints.
  const std::vector< double > hexa20Coords = {
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, // bottom corners
    0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, // top corners
    0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5, 0.0, // bottom-face edge midsides
    0.5, 0.0, 1.0, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0, 0.0, 0.5, 1.0, // top-face edge midsides
    0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5  // vertical edge midsides
  };
  const Eigen::Map< const Eigen::VectorXd > hexa20ParentCoords( hexa20Coords.data(), 60 );

  for ( int face = 1; face <= 6; face++ ) {
    BoundaryElement       boundaryEl( Hexa20, face, 3, hexa20ParentCoords );
    const Eigen::VectorXd Pk = boundaryEl.computeScalarLoadVector();

    Eigen::VectorXd parentScalar = Eigen::VectorXd::Zero( 20 );
    boundaryEl.assembleIntoParentScalar( Pk, parentScalar );

    throwExceptionOnFailure( checkIfEqual( parentScalar.sum(), 1.0, 1e-10 ),
                             "Hexa20 boundary computeScalarLoadVector() does not integrate to the face area "
                             "(face " +
                               std::to_string( face ) + ")." );
  }
}

// ---------------------------------------------------------------------------------------------
// Not-yet-implemented derivative methods
// ---------------------------------------------------------------------------------------------

void testComputeDScalarLoadVectorThrowsNotImplemented()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );

  bool threw = false;
  try {
    boundaryEl.computeDScalarLoadVector_dCoordinates();
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeDScalarLoadVector_dCoordinates() must throw (not yet implemented)." );
}

void testComputeDVectorialLoadVectorThrowsNotImplemented()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );

  bool threw = false;
  try {
    boundaryEl.computeDVectorialLoadVector_dCoordinates( Eigen::Vector2d( 1.0, 0.0 ) );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeDVectorialLoadVector_dCoordinates() must throw (not yet implemented)." );
}

// ---------------------------------------------------------------------------------------------
// computeDSurfaceNormalVectorialLoadVector_dCoordinates(): smoke test only. The analytic
// derivative is with respect to the boundary element's OWN (private, condensed) coordinate
// vector, which cannot be perturbed independently from outside without knowing the parent<->
// boundary index map, so a full numerical-differentiation cross-check is not attempted here --
// this only confirms the 2D and 3D branches run and return a finite, correctly-sized matrix.
// ---------------------------------------------------------------------------------------------

void testComputeDSurfaceNormalVectorialLoadVectorDCoordinatesReturnsFiniteMatrix()
{
  {
    const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
    BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );
    const Eigen::MatrixXd                     K = boundaryEl.computeDSurfaceNormalVectorialLoadVector_dCoordinates();
    throwExceptionOnFailure( K.rows() == 4 && K.cols() == 4,
                             "computeDSurfaceNormalVectorialLoadVector_dCoordinates() (2D) has the wrong size." );
    throwExceptionOnFailure( K.allFinite(),
                             "computeDSurfaceNormalVectorialLoadVector_dCoordinates() (2D) returned a "
                             "non-finite value." );
  }
  {
    const Eigen::Map< const Eigen::VectorXd > parentCoords( unitCubeCoords.data(), 24 );
    BoundaryElement                           boundaryEl( Hexa8, 1, 3, parentCoords );
    const Eigen::MatrixXd                     K = boundaryEl.computeDSurfaceNormalVectorialLoadVector_dCoordinates();
    throwExceptionOnFailure( K.rows() == 12 && K.cols() == 12,
                             "computeDSurfaceNormalVectorialLoadVector_dCoordinates() (3D) has the wrong size." );
    throwExceptionOnFailure( K.allFinite(),
                             "computeDSurfaceNormalVectorialLoadVector_dCoordinates() (3D) returned a "
                             "non-finite value." );
  }
}

// ---------------------------------------------------------------------------------------------
// condense*/assemble* round trips
// ---------------------------------------------------------------------------------------------

void testCondenseAndAssembleVectorialRoundTrip()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );

  const Eigen::VectorXd boundaryCoords = boundaryEl.condenseParentToBoundaryVectorial( parentCoords );

  Eigen::VectorXd reconstructed = Eigen::VectorXd::Zero( 8 );
  boundaryEl.assembleIntoParentVectorial( boundaryCoords, reconstructed );

  // Every reconstructed entry that is nonzero must match the original coordinate exactly, since
  // condense() extracted it verbatim from the same vector.
  for ( int i = 0; i < 8; i++ )
    if ( std::abs( reconstructed( i ) ) > 1e-12 )
      throwExceptionOnFailure( checkIfEqual( reconstructed( i ), parentCoords( i ), 1e-12 ),
                               "condenseParentToBoundaryVectorial()/assembleIntoParentVectorial() round trip "
                               "does not reproduce the original coordinate." );
}

void testCondenseAndAssembleScalarRoundTrip()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );

  // An arbitrary per-node scalar field on the 4 parent nodes.
  Eigen::VectorXd parentScalarField( 4 );
  parentScalarField << 10.0, 20.0, 30.0, 40.0;

  const Eigen::VectorXd boundaryScalar = boundaryEl.condenseParentToBoundaryScalar( parentScalarField );

  Eigen::VectorXd reconstructed = Eigen::VectorXd::Zero( 4 );
  boundaryEl.assembleIntoParentScalar( boundaryScalar, reconstructed );

  for ( int i = 0; i < 4; i++ )
    if ( std::abs( reconstructed( i ) ) > 1e-12 )
      throwExceptionOnFailure( checkIfEqual( reconstructed( i ), parentScalarField( i ), 1e-12 ),
                               "condenseParentToBoundaryScalar()/assembleIntoParentScalar() round trip does "
                               "not reproduce the original scalar value." );
}

// ---------------------------------------------------------------------------------------------
// assembleIntoParentStiffness{Scalar,Vectorial}(): the assembly negates the boundary matrix (see
// the "mind the negative sign" comment in the implementation).
// ---------------------------------------------------------------------------------------------

void testAssembleIntoParentStiffnessVectorialAppliesNegativeSign()
{
  const Eigen::Map< const Eigen::VectorXd > parentCoords( unitSquareCoords.data(), 8 );
  BoundaryElement                           boundaryEl( Quad4, 1, 2, parentCoords );

  const Eigen::VectorXd boundaryCoords = boundaryEl.condenseParentToBoundaryVectorial( parentCoords );
  const int             nBoundaryDof   = static_cast< int >( boundaryCoords.size() );

  const Eigen::MatrixXd KBoundary = Eigen::MatrixXd::Identity( nBoundaryDof, nBoundaryDof );
  Eigen::MatrixXd       KParent   = Eigen::MatrixXd::Zero( 8, 8 );
  boundaryEl.assembleIntoParentStiffnessVectorial( KBoundary, KParent );

  // The trace over the (unknown) mapped block must equal -trace(KBoundary), regardless of the
  // exact index mapping, since assembly only ever touches diagonal-preserving (i,i)->(map(i),map(i))
  // entries for an identity input.
  throwExceptionOnFailure( checkIfEqual( KParent.trace(), -KBoundary.trace(), 1e-12 ),
                           "assembleIntoParentStiffnessVectorial() does not apply the documented negative "
                           "sign." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testConstructorThrowsForUnsupportedParentShape,
    testComputeScalarLoadVectorIntegratesToFaceLengthQuad4,
    testComputeScalarLoadVectorIntegratesToFaceAreaHexa8,
    testComputeVectorialLoadVectorIntegratesToDirectionTimesLengthQuad4,
    testComputeVectorialLoadVectorIntegratesToDirectionTimesAreaHexa8,
    testComputeSurfaceNormalVectorialLoadVectorPointsOutwardQuad4,
    testComputeSurfaceNormalVectorialLoadVectorPointsOutwardHexa8,
    testComputeSurfaceNormalVectorialLoadVectorForAxisymmetricElementsScalesByTwoPiR,
    testQuad8AndHexa20BoundariesIntegrateToFaceLengthAndArea,
    testComputeDScalarLoadVectorThrowsNotImplemented,
    testComputeDVectorialLoadVectorThrowsNotImplemented,
    testComputeDSurfaceNormalVectorialLoadVectorDCoordinatesReturnsFiniteMatrix,
    testCondenseAndAssembleVectorialRoundTrip,
    testCondenseAndAssembleScalarRoundTrip,
    testAssembleIntoParentStiffnessVectorialAppliesNegativeSign,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
