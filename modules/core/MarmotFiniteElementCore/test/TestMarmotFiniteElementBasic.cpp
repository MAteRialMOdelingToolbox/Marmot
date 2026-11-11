#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotMassLumping.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <tuple>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Discretizations::Elements;

// ---------------------------------------------------------------------------------------------
// getElementShapeByMetric()
// ---------------------------------------------------------------------------------------------

void testGetElementShapeByMetricAllValidCombinations()
{
  throwExceptionOnFailure( getElementShapeByMetric( 1, 2 ) == Bar2, "nDim=1, nNodes=2 must be Bar2." );
  throwExceptionOnFailure( getElementShapeByMetric( 1, 3 ) == Bar3, "nDim=1, nNodes=3 must be Bar3." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 3 ) == Tri3, "nDim=2, nNodes=3 must be Tri3." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 4 ) == Quad4, "nDim=2, nNodes=4 must be Quad4." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 8 ) == Quad8, "nDim=2, nNodes=8 must be Quad8." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 9 ) == Quad9, "nDim=2, nNodes=9 must be Quad9." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 16 ) == Quad16, "nDim=2, nNodes=16 must be Quad16." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 4 ) == Tetra4, "nDim=3, nNodes=4 must be Tetra4." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 10 ) == Tetra10, "nDim=3, nNodes=10 must be Tetra10." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 8 ) == Hexa8, "nDim=3, nNodes=8 must be Hexa8." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 20 ) == Hexa20, "nDim=3, nNodes=20 must be Hexa20." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 27 ) == Hexa27, "nDim=3, nNodes=27 must be Hexa27." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 64 ) == Hexa64, "nDim=3, nNodes=64 must be Hexa64." );
}

void testGetElementShapeByMetricThrowsForInvalidCombinations()
{
  const std::vector< std::pair< int, int > > invalidCombinations = {
    { 1, 4 }, // no 4-node 1D shape
    { 2, 5 }, // no 5-node 2D shape
    { 3, 5 }, // no 5-node 3D shape
    { 4, 8 }, // no 4-D shapes at all
  };

  for ( const auto& [nDim, nNodes] : invalidCombinations ) {
    bool threw = false;
    try {
      getElementShapeByMetric( nDim, nNodes );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw,
                             "getElementShapeByMetric() must throw for nDim=" + std::to_string( nDim ) +
                               ", nNodes=" + std::to_string( nNodes ) + "." );
  }
}

// ---------------------------------------------------------------------------------------------
// NB() / Jacobian(): dynamic (runtime-sized) versions have no callers anywhere in the codebase --
// every element uses the templated (compile-time-sized) overloads declared alongside them in
// MarmotFiniteElement.h instead. Cross-validated against those already-trusted templated versions
// for a known Quad4 configuration, rather than re-deriving expected values independently.
// ---------------------------------------------------------------------------------------------

void testDynamicNBMatchesTemplatedNBForQuad4()
{
  Marmot::Discretizations::Elements::GeometryElement< 2, 4 > geo;
  const auto N = geo.N( Eigen::Vector2d( 0.2, 0.3 ) ); // NSized = Matrix<double,1,4>

  const auto templatedResult = Marmot::Discretizations::Elements::NB< 2, 4 >( N );
  const auto dynamicResult   = Marmot::Discretizations::Elements::NB( Eigen::VectorXd( N.transpose() ), 2 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( dynamicResult ), Eigen::MatrixXd( templatedResult ), 1e-12 ),
                           "The dynamic NB() must match the templated NB<nDim,nNodes>() for the same input." );
}

void testDynamicJacobianMatchesTemplatedJacobianForQuad4()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 2.0, 0.0, 3.0, 2.0, 0.0, 1.0 }; // general quadrilateral

  Marmot::Discretizations::Elements::GeometryElement< 2, 4 > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const Eigen::Vector2d xi( 0.2, 0.3 );
  const auto            dNdXi = geo.dNdXi( xi );

  const auto templatedJ = geo.Jacobian( dNdXi ); // uses the templated Jacobian<2,4> internally

  const Eigen::MatrixXd dynamicDNdXi = dNdXi;
  const Eigen::VectorXd dynamicCoords( Eigen::Map< const Eigen::VectorXd >( nodeCoordsVec.data(), 8 ) );
  const Eigen::MatrixXd dynamicJ = Marmot::Discretizations::Elements::Jacobian( dynamicDNdXi, dynamicCoords );

  throwExceptionOnFailure( checkIfEqual( dynamicJ, Eigen::MatrixXd( templatedJ ), 1e-12 ),
                           "The dynamic Jacobian() must match the templated Jacobian<nDim,nNodes>() for the same "
                           "input." );
}

// ---------------------------------------------------------------------------------------------
// getReferenceDomain() / getInterpolationOrder()
// ---------------------------------------------------------------------------------------------

void testReferenceDomainAndInterpolationOrderOfAllShapes()
{
  const std::vector< std::tuple< ElementShapes, ReferenceDomain, int > > expected = {
    { Bar2, ReferenceDomain::Line, 1 },
    { Bar3, ReferenceDomain::Line, 2 },
    { Tri3, ReferenceDomain::Triangle, 1 },
    { Quad4, ReferenceDomain::Quadrilateral, 1 },
    { Quad8, ReferenceDomain::Quadrilateral, 2 },
    { Quad9, ReferenceDomain::Quadrilateral, 2 },
    { Quad16, ReferenceDomain::Quadrilateral, 3 },
    { Tetra4, ReferenceDomain::Tetrahedron, 1 },
    { Tetra10, ReferenceDomain::Tetrahedron, 2 },
    { Hexa8, ReferenceDomain::Hexahedron, 1 },
    { Hexa20, ReferenceDomain::Hexahedron, 2 },
    { Hexa27, ReferenceDomain::Hexahedron, 2 },
    { Hexa64, ReferenceDomain::Hexahedron, 3 },
  };
  for ( const auto& [shape, domain, order] : expected ) {
    throwExceptionOnFailure( getReferenceDomain( shape ) == domain,
                             "getReferenceDomain() is wrong for shape " + std::to_string( shape ) + "." );
    throwExceptionOnFailure( getInterpolationOrder( shape ) == order,
                             "getInterpolationOrder() is wrong for shape " + std::to_string( shape ) + "." );
  }
}

void testNumberOfCornerNodesOfAllReferenceDomains()
{
  const std::vector< std::pair< ReferenceDomain, int > > expected = { { ReferenceDomain::Line, 2 },
                                                                      { ReferenceDomain::Triangle, 3 },
                                                                      { ReferenceDomain::Quadrilateral, 4 },
                                                                      { ReferenceDomain::Tetrahedron, 4 },
                                                                      { ReferenceDomain::Hexahedron, 8 } };
  for ( const auto& [domain, nCorners] : expected )
    throwExceptionOnFailure( getNumberOfCornerNodes( domain ) == nCorners,
                             "getNumberOfCornerNodes() is wrong for reference domain " +
                               std::to_string( static_cast< int >( domain ) ) + "." );
}

// ---------------------------------------------------------------------------------------------
// MassLumping::nLinearNodes(): the linear element blended into the manifold-based lumping has the corner nodes of
// the element; checked at compile time, since the callers use it as a template argument
// ---------------------------------------------------------------------------------------------

namespace {
  using Marmot::Discretizations::Elements::MassLumping::nLinearNodes;

  static_assert( nLinearNodes< 1, 2 >() == 2, "Bar2 is its own linear element." );
  static_assert( nLinearNodes< 1, 3 >() == 2, "The linear element of Bar3 is Bar2." );
  static_assert( nLinearNodes< 2, 3 >() == 3, "Tri3 is its own linear element." );
  static_assert( nLinearNodes< 2, 4 >() == 4, "Quad4 is its own linear element." );
  static_assert( nLinearNodes< 2, 8 >() == 4, "The linear element of Quad8 is Quad4." );
  static_assert( nLinearNodes< 2, 9 >() == 4, "The linear element of Quad9 is Quad4." );
  static_assert( nLinearNodes< 2, 16 >() == 4, "The linear element of Quad16 is Quad4." );
  static_assert( nLinearNodes< 3, 4 >() == 4, "Tetra4 is its own linear element." );
  static_assert( nLinearNodes< 3, 10 >() == 4, "The linear element of Tetra10 is Tetra4, not Hexa8." );
  static_assert( nLinearNodes< 3, 8 >() == 8, "Hexa8 is its own linear element." );
  static_assert( nLinearNodes< 3, 20 >() == 8, "The linear element of Hexa20 is Hexa8." );
  static_assert( nLinearNodes< 3, 27 >() == 8, "The linear element of Hexa27 is Hexa8." );
  static_assert( nLinearNodes< 3, 64 >() == 8, "The linear element of Hexa64 is Hexa8." );
} // namespace

// ---------------------------------------------------------------------------------------------
// Quadrature::getGaussPointInfo() / getNumGaussPoints(): the branches actually reached by every
// existing element and boundary test all use FullIntegration for Bar3/Tetra4/Tetra10 (or, for
// Bar2/Quad4, never use ReducedIntegration) and never construct a Tetra10 element (there is no
// registered element type using that shape) or pass an unsupported shape at all.
// ---------------------------------------------------------------------------------------------

void testGetGaussPointInfoReducedIntegrationBranches()
{
  using namespace Marmot::Discretizations::Elements::Quadrature;

  throwExceptionOnFailure( getGaussPointInfo( Bar2, ReducedIntegration ).size() == 1,
                           "Bar2 reduced integration must use the 1-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Bar3, ReducedIntegration ).size() == 2,
                           "Bar3 reduced integration must use the 2-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Quad4, ReducedIntegration ).size() == 1,
                           "Quad4 reduced integration must use the 1x1-point rule." );
}

void testGetGaussPointInfoTetrahedralShapes()
{
  using namespace Marmot::Discretizations::Elements::Quadrature;

  // Tetra4/Tetra10 quadrature rules don't distinguish integration type; both calls must succeed
  // and return a non-empty rule.
  throwExceptionOnFailure( getGaussPointInfo( Tetra4, FullIntegration ).size() > 0,
                           "Tetra4 must return a non-empty quadrature rule." );
  throwExceptionOnFailure( getGaussPointInfo( Tetra10, ReducedIntegration ).size() > 0,
                           "Tetra10 must return a non-empty quadrature rule." );
}

void testGetGaussPointInfoThrowsForUnsupportedShape()
{
  using namespace Marmot::Discretizations::Elements::Quadrature;

  bool threw = false;
  try {
    getGaussPointInfo( Quad9, FullIntegration );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "getGaussPointInfo() must throw for an unsupported shape (Quad9)." );
}

void testGetNumGaussPointsMatchesGaussPointInfoSize()
{
  using namespace Marmot::Discretizations::Elements::Quadrature;

  throwExceptionOnFailure( getNumGaussPoints( Hexa8, FullIntegration ) ==
                             static_cast< int >( getGaussPointInfo( Hexa8, FullIntegration ).size() ),
                           "getNumGaussPoints() must match getGaussPointInfo().size()." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testGetElementShapeByMetricAllValidCombinations,
    testGetElementShapeByMetricThrowsForInvalidCombinations,
    testDynamicNBMatchesTemplatedNBForQuad4,
    testDynamicJacobianMatchesTemplatedJacobianForQuad4,
    testReferenceDomainAndInterpolationOrderOfAllShapes,
    testNumberOfCornerNodesOfAllReferenceDomains,
    testGetGaussPointInfoReducedIntegrationBranches,
    testGetGaussPointInfoTetrahedralShapes,
    testGetGaussPointInfoThrowsForUnsupportedShape,
    testGetNumGaussPointsMatchesGaussPointInfoSize,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
