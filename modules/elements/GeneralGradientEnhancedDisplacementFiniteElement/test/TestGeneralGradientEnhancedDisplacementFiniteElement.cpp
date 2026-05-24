#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

void testBasicPropertiesQuad4PlaneStress()
{
  // Test static structural properties of the GCPS4 element:
  // 2D, 4-node quad, 1 nonlocal variable, plane stress
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 1;
  const int     elId          = 1;
  const auto    intType       = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( elId, intType, ElemType::SectionType::PlaneStress );

  throwExceptionOnFailure( element->getNNodes() == nNodes,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nNodes" );
  throwExceptionOnFailure( element->getNSpatialDimensions() == nDim,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDim" );

  // sizeLoadVector = nNodes*nDim + nNodes*nNonlocalVars = 4*2 + 4*1 = 12
  throwExceptionOnFailure( element->getNDofPerElement() == 12,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDof" );

  throwExceptionOnFailure( element->getElementShape() == "quad4",
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect element shape" );

  // Full integration of Quad4 uses 2x2 = 4 Gauss points
  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 4,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect number of QPs" );
}

void testBasicPropertiesQuad8PlaneStrain()
{
  // Test static structural properties of the GCPE8 element:
  // 2D, 8-node serendipity quad, 1 nonlocal variable, plane strain
  constexpr int nDim          = 2;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  const int     elId          = 2;
  const auto    intType       = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( elId, intType, ElemType::SectionType::PlaneStrain );

  throwExceptionOnFailure( element->getNNodes() == nNodes,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nNodes" );
  throwExceptionOnFailure( element->getNSpatialDimensions() == nDim,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDim" );

  // sizeLoadVector = 8*2 + 8*1 = 24
  throwExceptionOnFailure( element->getNDofPerElement() == 24,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDof" );

  throwExceptionOnFailure( element->getElementShape() == "quad8",
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect element shape" );

  // Full integration of Quad8 uses 3x3 = 9 Gauss points
  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 9,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect number of QPs" );
}

void testBasicPropertiesHex8Solid()
{
  // Test static structural properties of the GC3D8 element:
  // 3D, 8-node hex, 1 nonlocal variable, solid
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  const int     elId          = 3;
  const auto    intType       = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( elId, intType, ElemType::SectionType::Solid );

  throwExceptionOnFailure( element->getNNodes() == nNodes,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nNodes" );
  throwExceptionOnFailure( element->getNSpatialDimensions() == nDim,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDim" );

  // sizeLoadVector = 8*3 + 8*1 = 32
  throwExceptionOnFailure( element->getNDofPerElement() == 32,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect nDof" );

  throwExceptionOnFailure( element->getElementShape() == "hexa8",
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect element shape" );

  // Full integration of Hex8 uses 2x2x2 = 8 Gauss points
  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 8,
                           MakeString() << __PRETTY_FUNCTION__ << ": incorrect number of QPs" );
}

void testDofIndicesPermutationPatternQuad4()
{
  // Permutation pattern for GCPS4: maps from split [u-block | k-block] to
  // interleaved [u1x, u1y, k1, u2x, u2y, k2, ...] DOF layout used in the solver
  // For nDim=2, nNodes=4, nNonlocalVars=1, nNonLocalNodes=4:
  // Expected: {0, 1, 3, 4, 6, 7, 9, 10, 2, 5, 8, 11}
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 1;
  const int     elId          = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( elId,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStress );

  const std::vector< int > pattern  = element->getDofIndicesPermutationPattern();
  const std::vector< int > expected = { 0, 1, 3, 4, 6, 7, 9, 10, 2, 5, 8, 11 };

  throwExceptionOnFailure( pattern.size() == expected.size(),
                           MakeString() << __PRETTY_FUNCTION__ << ": pattern size mismatch" );
  for ( std::size_t i = 0; i < expected.size(); ++i ) {
    throwExceptionOnFailure( pattern[i] == expected[i],
                             MakeString() << __PRETTY_FUNCTION__ << ": pattern[" << i << "] mismatch" );
  }
}

void testDofIndicesPermutationPatternQuad4TwoNonlocalVars()
{
  // Permutation pattern for G2GCPS4: nDim=2, nNodes=4, nNonlocalVars=2
  // Expected: {0, 1, 4, 5, 8, 9, 12, 13, 2, 6, 10, 14, 3, 7, 11, 15}
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 2;
  const int     elId          = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( elId,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStress );

  const std::vector< int > pattern  = element->getDofIndicesPermutationPattern();
  const std::vector< int > expected = { 0, 1, 4, 5, 8, 9, 12, 13, 2, 6, 10, 14, 3, 7, 11, 15 };

  throwExceptionOnFailure( pattern.size() == expected.size(),
                           MakeString() << __PRETTY_FUNCTION__ << ": pattern size mismatch" );
  for ( std::size_t i = 0; i < expected.size(); ++i ) {
    throwExceptionOnFailure( pattern[i] == expected[i],
                             MakeString() << __PRETTY_FUNCTION__ << ": pattern[" << i << "] mismatch" );
  }
}

void testNodeFieldsOneNonlocalVar()
{
  // Node fields for nDim=2, nNodes=4, nNonlocalVars=1:
  // every node should have {"displacement", "nonlocal damage"}
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStress );

  const auto nodeFields = element->getNodeFields();
  throwExceptionOnFailure( static_cast< int >( nodeFields.size() ) == nNodes,
                           MakeString() << __PRETTY_FUNCTION__ << ": nodeFields size mismatch" );
  for ( int i = 0; i < nNodes; ++i ) {
    throwExceptionOnFailure( nodeFields[i].size() == 2,
                             MakeString() << __PRETTY_FUNCTION__ << ": node " << i << " should have 2 fields" );
    throwExceptionOnFailure( nodeFields[i][0] == "displacement",
                             MakeString() << __PRETTY_FUNCTION__ << ": node " << i << " first field != displacement" );
    throwExceptionOnFailure( nodeFields[i][1] == "nonlocal damage",
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": node " << i << " second field != nonlocal damage" );
  }
}

void testCoordinatesAtCenter()
{
  // For a unit square element [0,1]x[0,1], the center should be at (0.5, 0.5)
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStress );

  const std::vector< double > nodeCoords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  element->assignNodeCoordinates( nodeCoords.data() );

  const auto centerCoords = element->getCoordinatesAtCenter();
  throwExceptionOnFailure( static_cast< int >( centerCoords.size() ) == nDim,
                           MakeString() << __PRETTY_FUNCTION__ << ": wrong dimension of center coordinates" );
  throwExceptionOnFailure( checkIfEqual( centerCoords[0], 0.5 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": wrong x-coordinate at center" );
  throwExceptionOnFailure( checkIfEqual( centerCoords[1], 0.5 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": wrong y-coordinate at center" );
}

void testCoordinatesAtQuadraturePoints()
{
  // For a unit square, all 4 QPs must lie strictly inside [0,1]x[0,1]
  constexpr int nDim          = 2;
  constexpr int nNodes        = 4;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStress );

  const std::vector< double > nodeCoords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  element->assignNodeCoordinates( nodeCoords.data() );

  const auto qpCoords = element->getCoordinatesAtQuadraturePoints();
  throwExceptionOnFailure( static_cast< int >( qpCoords.size() ) == element->getNumberOfQuadraturePoints(),
                           MakeString() << __PRETTY_FUNCTION__ << ": wrong number of QP coordinate sets" );

  for ( const auto& coords : qpCoords ) {
    throwExceptionOnFailure( static_cast< int >( coords.size() ) == nDim,
                             MakeString() << __PRETTY_FUNCTION__ << ": wrong QP coordinate dimension" );
    throwExceptionOnFailure( coords[0] > 0.0 && coords[0] < 1.0,
                             MakeString() << __PRETTY_FUNCTION__ << ": QP x-coordinate out of element bounds" );
    throwExceptionOnFailure( coords[1] > 0.0 && coords[1] < 1.0,
                             MakeString() << __PRETTY_FUNCTION__ << ": QP y-coordinate out of element bounds" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testBasicPropertiesQuad4PlaneStress,
    testBasicPropertiesQuad8PlaneStrain,
    testBasicPropertiesHex8Solid,
    testDofIndicesPermutationPatternQuad4,
    testDofIndicesPermutationPatternQuad4TwoNonlocalVars,
    testNodeFieldsOneNonlocalVar,
    testCoordinatesAtCenter,
    testCoordinatesAtQuadraturePoints,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
