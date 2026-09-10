#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <algorithm>
#include <string>

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

// ---------------------------------------------------------------------------------------------
// Lumped (diagonal) mass/capacity matrix tests
//
// computeLumpedInertia() and computeLumpedDamping() both use the manifold-based scheme of Yang
// et al. (2017) for their non-zero block -- computeLumpedInertia() for the displacement block
// (using density), computeLumpedDamping() for the non-local field block (using non-local
// viscosity), mixing the high-order shape function N with the corresponding corner-node linear
// shape function N_lin via N_weighted = w*N + (1-w)*N_lin (only the corner entries receive the
// N_lin correction), with w = 1/2 by default and w = 1/3 special-cased for Hexa20: with the
// default split, the negative corner contribution of the Hexa20 serendipity shape function
// exactly cancels the positive corner contribution of the trilinear shape function for any
// regular element, producing exactly zero corner mass (see TestDisplacementFiniteElement.cpp for
// the derivation).
//
// Both tests below use equal-order interpolation (nNonLocalNodes == nNodes, the default), with
// density == non-local viscosity == 1, so the non-local block reduces to the exact same
// shape-function diagonal pattern as the displacement block, and the reference values (from an
// independent SymPy computation) apply identically to both blocks. No micro-inertia is assigned,
// so computeLumpedInertia()'s non-local block is exactly zero -- checked here too.
// ---------------------------------------------------------------------------------------------

void checkQuad8AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes intType, const std::string& label )
{
  constexpr int nDim          = 2;
  constexpr int nNodes        = 8; // Quad8 (quadratic serendipity)
  constexpr int nNonlocalVars = 1;
  const int     elId          = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;
  const auto secType          = ElemType::SectionType::PlaneStrain;

  // Unit square with midside nodes at exact edge midpoints (straight edges). Node ordering
  // per MarmotFiniteElement2D.cpp Quad8::N: 0-3 corners CCW, 4-7 midsides.
  const std::vector< double > nodeCoordsVec = { 0.0,
                                                0.0,
                                                1.0,
                                                0.0,
                                                1.0,
                                                1.0,
                                                0.0,
                                                1.0, // corners
                                                0.5,
                                                0.0,
                                                1.0,
                                                0.5,
                                                0.5,
                                                1.0,
                                                0.0,
                                                0.5 }; // midsides

  auto element = std::make_unique< ElemType >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  // AT2PhaseField properties: E, nu, Gc, l, density, nonlocalViscosity
  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 }; // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  // Reference values (see TestDisplacementFiniteElement.cpp for the SymPy derivation): corners
  // = 1/12, midsides = 1/6. Because the geometry map is affine here, both full (3x3) and reduced
  // (2x2) Gauss integrate this exactly, so both integration types reproduce the same values.
  const double expectedCorner  = 1.0 / 12.0;
  const double expectedMidside = 1.0 / 6.0;

  constexpr int sizeDoFU = nNodes * nDim;

  double totalMass = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 4 ? expectedCorner : expectedMidside;
    throwExceptionOnFailure( M[i * nDim] > 0.0,
                             label + ": Quad8 displacement lumped mass entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], expected, 1e-10 ),
                             label +
                               ": Quad8 displacement lumped mass entry does not match the analytic reference "
                               "(node " +
                               std::to_string( i ) + ")." );
    totalMass += M[i * nDim];
  }
  throwExceptionOnFailure( checkIfEqual( totalMass, 1.0, 1e-10 ),
                           label + ": Quad8 displacement block does not conserve the total element mass." );

  // No micro-inertia was assigned, so computeLumpedInertia()'s non-local block must be exactly
  // zero -- the field stays first order in time.
  for ( int i = 0; i < nNodes; i++ )
    throwExceptionOnFailure( M[sizeDoFU + i] == 0.0,
                             label +
                               ": Quad8 non-local lumped inertia entry is non-zero without a micro-inertia "
                               "assigned (node " +
                               std::to_string( i ) + ")." );

  std::vector< double > C( element->getNDofPerElement(), 0.0 );
  element->computeLumpedDamping( C.data() );

  double totalCapacity = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 4 ? expectedCorner : expectedMidside;
    throwExceptionOnFailure( C[sizeDoFU + i] > 0.0,
                             label + ": Quad8 non-local lumped damping entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    throwExceptionOnFailure( checkIfEqual( C[sizeDoFU + i], expected, 1e-10 ),
                             label +
                               ": Quad8 non-local lumped damping entry does not match the analytic "
                               "reference (node " +
                               std::to_string( i ) + ")." );
    totalCapacity += C[sizeDoFU + i];
  }
  throwExceptionOnFailure( checkIfEqual( totalCapacity, 1.0, 1e-10 ),
                           label + ": Quad8 non-local block does not conserve the total element capacity." );
}

void testLumpedInertiaQuad8FullIntegrationMatchesAnalyticValues()
{
  checkQuad8AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes::FullIntegration, "FullIntegration" );
}

void testLumpedInertiaQuad8ReducedIntegrationMatchesAnalyticValues()
{
  checkQuad8AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                                  "ReducedIntegration" );
}

void checkHexa20AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes intType, const std::string& label )
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 20; // Hexa20 (quadratic serendipity)
  constexpr int nNonlocalVars = 1;
  const int     elId          = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;
  const auto secType          = ElemType::SectionType::Solid;

  // Unit cube with edge-midside nodes at exact midpoints (straight edges). Node ordering per
  // MarmotFiniteElement3D.cpp Hexa20::N: 0-7 corners, 8-19 edge midsides.
  const std::vector< double > nodeCoordsVec = { // corners
                                                0.0,
                                                0.0,
                                                0.0,
                                                1.0,
                                                0.0,
                                                0.0,
                                                1.0,
                                                1.0,
                                                0.0,
                                                0.0,
                                                1.0,
                                                0.0,
                                                0.0,
                                                0.0,
                                                1.0,
                                                1.0,
                                                0.0,
                                                1.0,
                                                1.0,
                                                1.0,
                                                1.0,
                                                0.0,
                                                1.0,
                                                1.0,
                                                // bottom-face edge midsides (0-1, 1-2, 2-3, 3-0)
                                                0.5,
                                                0.0,
                                                0.0,
                                                1.0,
                                                0.5,
                                                0.0,
                                                0.5,
                                                1.0,
                                                0.0,
                                                0.0,
                                                0.5,
                                                0.0,
                                                // top-face edge midsides (4-5, 5-6, 6-7, 7-4)
                                                0.5,
                                                0.0,
                                                1.0,
                                                1.0,
                                                0.5,
                                                1.0,
                                                0.5,
                                                1.0,
                                                1.0,
                                                0.0,
                                                0.5,
                                                1.0,
                                                // vertical edge midsides (0-4, 1-5, 2-6, 3-7)
                                                0.0,
                                                0.0,
                                                0.5,
                                                1.0,
                                                0.0,
                                                0.5,
                                                1.0,
                                                1.0,
                                                0.5,
                                                0.0,
                                                1.0,
                                                0.5 };

  auto element = std::make_unique< ElemType >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  // AT2PhaseField properties: E, nu, Gc, l, density, nonlocalViscosity
  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );

  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  // Reference values (see TestDisplacementFiniteElement.cpp for the SymPy derivation of the
  // Hexa20 1/3-2/3 special case): corners = 1/24, edges = 1/18.
  const double expectedCorner = 1.0 / 24.0;
  const double expectedEdge   = 1.0 / 18.0;

  constexpr int sizeDoFU = nNodes * nDim;

  double totalMass = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 8 ? expectedCorner : expectedEdge;
    throwExceptionOnFailure( M[i * nDim] > 0.0,
                             label + ": Hexa20 displacement lumped mass entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], expected, 1e-10 ),
                             label +
                               ": Hexa20 displacement lumped mass entry does not match the analytic "
                               "reference (node " +
                               std::to_string( i ) + ")." );
    totalMass += M[i * nDim];
  }
  throwExceptionOnFailure( checkIfEqual( totalMass, 1.0, 1e-10 ),
                           label + ": Hexa20 displacement block does not conserve the total element mass." );

  // No micro-inertia was assigned, so computeLumpedInertia()'s non-local block must be exactly
  // zero -- the field stays first order in time.
  for ( int i = 0; i < nNodes; i++ )
    throwExceptionOnFailure( M[sizeDoFU + i] == 0.0,
                             label +
                               ": Hexa20 non-local lumped inertia entry is non-zero without a micro-inertia "
                               "assigned (node " +
                               std::to_string( i ) + ")." );

  std::vector< double > C( element->getNDofPerElement(), 0.0 );
  element->computeLumpedDamping( C.data() );

  double totalCapacity = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 8 ? expectedCorner : expectedEdge;
    throwExceptionOnFailure( C[sizeDoFU + i] > 0.0,
                             label + ": Hexa20 non-local lumped damping entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    throwExceptionOnFailure( checkIfEqual( C[sizeDoFU + i], expected, 1e-10 ),
                             label +
                               ": Hexa20 non-local lumped damping entry does not match the analytic "
                               "reference (node " +
                               std::to_string( i ) + ")." );
    totalCapacity += C[sizeDoFU + i];
  }
  throwExceptionOnFailure( checkIfEqual( totalCapacity, 1.0, 1e-10 ),
                           label + ": Hexa20 non-local block does not conserve the total element capacity." );
}

void testLumpedInertiaHexa20FullIntegrationMatchesAnalyticValues()
{
  checkHexa20AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes::FullIntegration, "FullIntegration" );
}

void testLumpedInertiaHexa20ReducedIntegrationMatchesAnalyticValues()
{
  checkHexa20AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                                   "ReducedIntegration" );
}

/**
 * @brief Build a cube of the given edge length with AT2PhaseField and, optionally, a micro-inertia.
 * @param edgeLength Edge length of the cube.
 * @param microInertia Micro-inertia to assign, in seconds squared; not assigned at all if zero.
 * @param nonlocalViscosity The material's non-local viscosity, which is the damping of the
 *        non-local field once that field carries a micro-inertia.
 * @param density Material density; raising it slows the mechanical waves and so pushes the
 *        mechanical limit out of the way of the non-local one.
 * @return The critical time step the element reports.
 */
double criticalTimeStepOfCube( double edgeLength, double microInertia, double nonlocalViscosity, double density )
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  const int     elId          = 1;
  const auto    intType       = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  const double                h             = edgeLength;
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, h, 0.0, 0.0, h, h, 0.0, 0.0, h, 0.0,
                                                0.0, 0.0, h,   h, 0.0, h,   h, h, h,   0.0, h, h };

  auto element = std::make_unique< ElemType >( elId, intType, ElemType::SectionType::Solid );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  // AT2PhaseField properties: E, nu, Gc, l, density, nonlocalViscosity. The internal length l
  // enters the element as c = l^2, and l >> h is what puts the non-local field in the regime where
  // its limit is set by the Laplacian rather than by the reaction term.
  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, density, nonlocalViscosity };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  if ( microInertia > 0.0 )
    element->assignProperty( "nonlocal micro inertia", &microInertia, 1 );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const std::vector< double > QTotal( element->getNDofPerElement(), 0.0 );

  double criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  return criticalTimeStep;
}

void testNonlocalMicroInertiaIsZeroUnlessAssigned()
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::Solid );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  // Deliberately dirtied: the contract is that the element WRITES the whole vector, so an element
  // without a micro-inertia has to zero its non-local block rather than leave what happened to be
  // there. The displacement block is not checked against zero here: it legitimately carries the
  // density-based mass computeLumpedInertia() now reports for every field at once.
  constexpr int         sizeDoFU = nNodes * nDim;
  std::vector< double > microInertia( element->getNDofPerElement(), 12345.0 );
  element->computeLumpedInertia( microInertia.data() );

  for ( size_t i = sizeDoFU; i < microInertia.size(); i++ )
    throwExceptionOnFailure( microInertia[i] == 0.0,
                             MakeString()
                               << __PRETTY_FUNCTION__
                               << ": an element without the property reported a micro-inertia at dof " << i );
}

void testNonlocalMicroInertiaSumsToMicroInertiaTimesVolume()
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  const double                h             = 2.0;
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, h, 0.0, 0.0, h, h, 0.0, 0.0, h, 0.0,
                                                0.0, 0.0, h,   h, 0.0, h,   h, h, h,   0.0, h, h };

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::Solid );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const double microInertiaValue = 3.0e-7;
  element->assignProperty( "nonlocal micro inertia", &microInertiaValue, 1 );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  // The displacement block carries the ordinary density-based mass now, not the micro-inertia --
  // that only ever lands in the non-local block, which is what is summed below.
  constexpr int         sizeDoFU = nNodes * nDim;
  std::vector< double > microInertia( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( microInertia.data() );

  double total = 0.0;
  for ( size_t i = sizeDoFU; i < microInertia.size(); i++ )
    total += microInertia[i];

  // Row-summing a consistent mass matrix conserves its total exactly, whatever the lumping weights
  // do between nodes, so this is an exact statement and not a tolerance-fitted one. It also fails
  // loudly on a mis-ordered node list, which would fold the cube inside out and change its volume.
  throwExceptionOnFailure( checkIfEqual( total, microInertiaValue * h * h * h, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": the lumped micro-inertia sums to " << total
                                        << " instead of " << microInertiaValue * h * h * h );
}

void testNonlocalMicroInertiaPropertyIsValidated()
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::Solid );

  const auto rejects = [&]( const double* values, int nValues ) {
    try {
      element->assignProperty( "nonlocal micro inertia", values, nValues );
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  };

  const std::vector< double > twoValues = { 1.0e-12, 1.0e-12 };
  throwExceptionOnFailure( rejects( twoValues.data(), 2 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << ": one value per non-local variable was not enforced" );

  const double negative = -1.0e-12;
  throwExceptionOnFailure( rejects( &negative, 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": a negative micro-inertia was accepted" );
}

void testBulkViscosityDamageDegradationIsValidated()
{
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::Solid );

  // reached through the base pointer, the way the EdelweissFE wrapper reaches it
  MarmotElement* base = element.get();

  const auto rejects = [&]( const double* values, int nValues ) {
    try {
      base->assignProperty( "bulk viscosity damage degradation", values, nValues );
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  };

  const std::vector< double > twoValues = { 2.0, 2.0 };
  throwExceptionOnFailure( rejects( twoValues.data(), 2 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": the single-value arity was not enforced" );

  const double negative = -1.0;
  throwExceptionOnFailure( rejects( &negative, 1 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << ": a negative exponent was accepted, which would amplify the viscous stress "
                                           "as "
                                           "the material fails" );

  const double exponent = 2.0;
  throwExceptionOnFailure( !rejects( &exponent, 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": a valid exponent was rejected" );

  // The element must also report the name, or a deck asking for it would be silently ignored by
  // whatever validates against getPropertyNames().
  const auto names = base->getPropertyNames();
  throwExceptionOnFailure( std::find( names.begin(), names.end(), "bulk viscosity damage degradation" ) != names.end(),
                           MakeString() << __PRETTY_FUNCTION__
                                        << ": the degradation is assignable but not reported by getPropertyNames()" );
}

void testNonlocalCriticalTimeStepScalesWithSqrtOfMicroInertia()
{
  // No viscosity, so the damping factor is exactly one and the limit is exactly 2 / omega_max --
  // which makes the expected ratio exact rather than approximate. A high density keeps the
  // mechanical limit far above the non-local one, so what is measured is the non-local branch.
  const double density = 1.0e4;

  const double dtBase      = criticalTimeStepOfCube( 0.2, 1.0e-4, 0.0, density );
  const double dtFourfold  = criticalTimeStepOfCube( 0.2, 4.0e-4, 0.0, density );
  const double dtMechanica = criticalTimeStepOfCube( 0.2, 0.0, 0.0, density );

  throwExceptionOnFailure( dtBase < dtMechanica,
                           MakeString() << __PRETTY_FUNCTION__
                                        << ": the non-local field did not govern, so the ratio below would be "
                                           "measuring the mechanical limit instead" );

  // A second-order field's limit goes with the square root of its inertia. A first-order (viscous)
  // one would go with the coefficient itself, so this ratio is what separates the two.
  throwExceptionOnFailure( checkIfEqual( dtFourfold / dtBase, 2.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": expected a ratio of 2, got "
                                        << dtFourfold / dtBase );
}

void testNonlocalCriticalTimeStepScalesLinearlyWithElementSize()
{
  const double density = 1.0e4;

  // h << l on both meshes, so the Laplacian dominates the reaction term and the limit is in its
  // asymptotic regime.
  const double dtCoarse = criticalTimeStepOfCube( 0.2, 1.0e-4, 0.0, density );
  const double dtFine   = criticalTimeStepOfCube( 0.1, 1.0e-4, 0.0, density );

  throwExceptionOnFailure( dtFine < criticalTimeStepOfCube( 0.1, 0.0, 0.0, density ),
                           MakeString() << __PRETTY_FUNCTION__ << ": the non-local field did not govern" );

  // This is the whole point of the micro-inertia: halving the element size halves the stable
  // increment. The parabolic scheme it replaces would have quartered it, so 0.25 is what this test
  // exists to exclude.
  const double ratio = dtFine / dtCoarse;
  throwExceptionOnFailure( std::abs( ratio - 0.5 ) < 0.01,
                           MakeString() << __PRETTY_FUNCTION__ << ": expected a ratio of 0.5, got " << ratio );
}

void testTinyMicroInertiaIsNotTreatedAsAbsent()
{
  /* A micro-inertia is a time SQUARED, so ordinary values are very small: eta^2/4 is 2.5e-13 for a
   * viscosity of 1e-6 s. Eigen's isZero() tests against a precision threshold of 1e-12 by default,
   * which would call such a value absent -- the property assigned, the field integrated with no
   * second-order term, and its stability limit unchecked, all silently. Hence an exact comparison,
   * and hence this test, whose value sits two orders BELOW that threshold.
   */
  const double tinyMicroInertia = 1.0e-14;
  const double density          = 1.0e4;

  const double dtWithTiny   = criticalTimeStepOfCube( 0.2, tinyMicroInertia, 0.0, density );
  const double dtMechanical = criticalTimeStepOfCube( 0.2, 0.0, 0.0, density );

  throwExceptionOnFailure( dtWithTiny < dtMechanical,
                           MakeString() << __PRETTY_FUNCTION__ << ": a micro-inertia of " << tinyMicroInertia
                                        << " was ignored -- the stable increment stayed at the mechanical "
                                        << dtMechanical );

  // And it reaches the assembled vector, not just the time step.
  constexpr int nDim          = 3;
  constexpr int nNodes        = 8;
  constexpr int nNonlocalVars = 1;
  using ElemType              = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVars >;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::Solid );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 20000.0, 0.2, 1.0, 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "AT2PHASEFIELD", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );
  element->assignProperty( materialSection );
  element->assignProperty( "nonlocal micro inertia", &tinyMicroInertia, 1 );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > microInertia( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( microInertia.data() );

  double total = 0.0;
  for ( size_t i = nNodes * nDim; i < microInertia.size(); i++ )
    total += microInertia[i];

  throwExceptionOnFailure( checkIfEqual( total, tinyMicroInertia, 1e-24 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": the lumped micro-inertia sums to " << total
                                        << " instead of " << tinyMicroInertia );
}

void testNonlocalViscosityLowersTheCriticalTimeStep()
{
  const double density = 1.0e4;

  const double dtUndamped = criticalTimeStepOfCube( 0.2, 1.0e-4, 0.0, density );
  const double dtDamped   = criticalTimeStepOfCube( 0.2, 1.0e-4, 1.0e-2, density );

  // Damping never buys time step under central differences, it costs it. Getting the sign of that
  // wrong would leave a damped run integrating above its limit while the estimate reported safety.
  throwExceptionOnFailure( dtDamped < dtUndamped,
                           MakeString() << __PRETTY_FUNCTION__ << ": damping did not lower the stable increment ("
                                        << dtDamped << " against " << dtUndamped << ")" );
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
    testLumpedInertiaQuad8FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaQuad8ReducedIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20ReducedIntegrationMatchesAnalyticValues,
    testNonlocalMicroInertiaIsZeroUnlessAssigned,
    testNonlocalMicroInertiaSumsToMicroInertiaTimesVolume,
    testNonlocalMicroInertiaPropertyIsValidated,
    testBulkViscosityDamageDegradationIsValidated,
    testNonlocalCriticalTimeStepScalesWithSqrtOfMicroInertia,
    testNonlocalCriticalTimeStepScalesLinearlyWithElementSize,
    testNonlocalViscosityLowersTheCriticalTimeStep,
    testTinyMicroInertiaIsNotTreatedAsAbsent,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
