#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement;

// ---------------------------------------------------------------------------------------------
// getElementShapeByMetric()
// ---------------------------------------------------------------------------------------------

void testGetElementShapeByMetricAllValidCombinations()
{
  throwExceptionOnFailure( getElementShapeByMetric( 1, 2 ) == Bar2, "nDim=1, nNodes=2 must be Bar2." );
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
    { 1, 3 }, // no 3-node 1D shape
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
  MarmotGeometryElement< 2, 4 > geo;
  const auto                    N = geo.N( Eigen::Vector2d( 0.2, 0.3 ) ); // NSized = Matrix<double,1,4>

  const auto templatedResult = Marmot::FiniteElement::NB< 2, 4 >( N );
  const auto dynamicResult   = Marmot::FiniteElement::NB( Eigen::VectorXd( N.transpose() ), 2 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( dynamicResult ), Eigen::MatrixXd( templatedResult ), 1e-12 ),
                           "The dynamic NB() must match the templated NB<nDim,nNodes>() for the same input." );
}

void testDynamicJacobianMatchesTemplatedJacobianForQuad4()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 2.0, 0.0, 3.0, 2.0, 0.0, 1.0 }; // general quadrilateral

  MarmotGeometryElement< 2, 4 > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const Eigen::Vector2d xi( 0.2, 0.3 );
  const auto            dNdXi = geo.dNdXi( xi );

  const auto templatedJ = geo.Jacobian( dNdXi ); // uses the templated Jacobian<2,4> internally

  const Eigen::MatrixXd dynamicDNdXi = dNdXi;
  const Eigen::VectorXd dynamicCoords( Eigen::Map< const Eigen::VectorXd >( nodeCoordsVec.data(), 8 ) );
  const Eigen::MatrixXd dynamicJ = Marmot::FiniteElement::Jacobian( dynamicDNdXi, dynamicCoords );

  throwExceptionOnFailure( checkIfEqual( dynamicJ, Eigen::MatrixXd( templatedJ ), 1e-12 ),
                           "The dynamic Jacobian() must match the templated Jacobian<nDim,nNodes>() for the same "
                           "input." );
}

// ---------------------------------------------------------------------------------------------
// Quadrature::getGaussPointInfo() / getNumGaussPoints(): the branches actually reached by every
// existing element and boundary test all use FullIntegration for Bar3/Tetra4/Tetra10 (or, for
// Bar2/Quad4, never use ReducedIntegration) and never construct a Tetra4/Tetra10 element (there is
// no registered element type using either shape) or pass an unsupported shape at all.
// ---------------------------------------------------------------------------------------------

void testGetGaussPointInfoReducedIntegrationBranches()
{
  using namespace Marmot::FiniteElement::Quadrature;

  throwExceptionOnFailure( getGaussPointInfo( Bar2, ReducedIntegration ).size() == 1,
                           "Bar2 reduced integration must use the 1-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Bar3, ReducedIntegration ).size() == 2,
                           "Bar3 reduced integration must use the 2-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Quad4, ReducedIntegration ).size() == 1,
                           "Quad4 reduced integration must use the 1x1-point rule." );
}

void testGetGaussPointInfoTetrahedralShapes()
{
  using namespace Marmot::FiniteElement::Quadrature;

  // Tetra4/Tetra10 quadrature rules don't distinguish integration type; both calls must succeed
  // and return a non-empty rule.
  throwExceptionOnFailure( getGaussPointInfo( Tetra4, FullIntegration ).size() > 0,
                           "Tetra4 must return a non-empty quadrature rule." );
  throwExceptionOnFailure( getGaussPointInfo( Tetra10, ReducedIntegration ).size() > 0,
                           "Tetra10 must return a non-empty quadrature rule." );
}

void testGetGaussPointInfoThrowsForUnsupportedShape()
{
  using namespace Marmot::FiniteElement::Quadrature;

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
  using namespace Marmot::FiniteElement::Quadrature;

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
    testGetGaussPointInfoReducedIntegrationBranches,
    testGetGaussPointInfoTetrahedralShapes,
    testGetGaussPointInfoThrowsForUnsupportedShape,
    testGetNumGaussPointsMatchesGaussPointInfoSize,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
