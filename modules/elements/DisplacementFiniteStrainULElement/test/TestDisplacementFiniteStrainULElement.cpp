#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

// ---------------------------------------------------------------------------------------------
// Lumped (diagonal) mass matrix tests
//
// computeLumpedInertia() uses the manifold-based scheme of Yang et al. (2017), mixing the
// high-order shape function N with the corresponding corner-node linear shape function N_lin
// via N_weighted = w*N + (1-w)*N_lin (only the corner entries receive the N_lin correction),
// with w = 1/2 by default and w = 1/3 special-cased for Hexa20: with the default split, the
// negative corner contribution of the Hexa20 serendipity shape function exactly cancels the
// positive corner contribution of the trilinear shape function for any regular element,
// producing exactly zero corner mass (see TestDisplacementFiniteElement.cpp for the derivation).
// The reference values below come from an independent SymPy computation; they are identical to
// the DisplacementFiniteElement case because both classes share the same Quad8/Hexa20 shape
// functions via MarmotGeometryElement<nDim,nNodes>.
// ---------------------------------------------------------------------------------------------

void testLumpedInertiaHexa8RegularElementIsPositiveAndConservesMass()
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 8; // Hexa8 (linear)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteStrainULElement< nDim, nNodes >::SectionType::Solid;

  // Unit cube; node ordering per MarmotFiniteElement3D.cpp Hexa8::N.
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 1.0, 1.0, density }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );

  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  // Unit cube, uniform density: by symmetry each corner carries exactly volume*density/8.
  for ( int i = 0; i < nNodes; i++ ) {
    throwExceptionOnFailure( M[i * nDim] > 0.0, "Hexa8 lumped mass entry is not strictly positive." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], 0.125, 1e-12 ),
                             "Hexa8 lumped mass entry does not match the analytic value." );
  }
}

void checkQuad8AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes intType, const std::string& label )
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 8; // Quad8 (quadratic serendipity)
  const int     elId    = 1;
  const auto    secType = DisplacementFiniteStrainULElement< nDim, nNodes >::SectionType::PlaneStrain;

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

  auto element = std::make_unique< DisplacementFiniteStrainULElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 1.0, 1.0, density }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };             // thickness
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

  double totalMass = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 4 ? expectedCorner : expectedMidside;
    throwExceptionOnFailure( M[i * nDim] > 0.0,
                             label + ": Quad8 lumped mass entry is not strictly positive (node " + std::to_string( i ) +
                               ")." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], expected, 1e-10 ),
                             label + ": Quad8 lumped mass entry does not match the analytic reference (node " +
                               std::to_string( i ) + ")." );
    totalMass += M[i * nDim];
  }
  throwExceptionOnFailure( checkIfEqual( totalMass, 1.0, 1e-10 ),
                           label + ": Quad8 lumped mass does not conserve the total element mass." );
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
  constexpr int nDim    = 3;
  constexpr int nNodes  = 20; // Hexa20 (quadratic serendipity)
  const int     elId    = 1;
  const auto    secType = DisplacementFiniteStrainULElement< nDim, nNodes >::SectionType::Solid;

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

  auto element = std::make_unique< DisplacementFiniteStrainULElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 1.0, 1.0, density }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );

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

  double totalMass = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    const double expected = i < 8 ? expectedCorner : expectedEdge;
    throwExceptionOnFailure( M[i * nDim] > 0.0,
                             label + ": Hexa20 lumped mass entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], expected, 1e-10 ),
                             label + ": Hexa20 lumped mass entry does not match the analytic reference (node " +
                               std::to_string( i ) + ")." );
    totalMass += M[i * nDim];
  }
  throwExceptionOnFailure( checkIfEqual( totalMass, 1.0, 1e-10 ),
                           label + ": Hexa20 lumped mass does not conserve the total element mass." );
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

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testLumpedInertiaHexa8RegularElementIsPositiveAndConservesMass,
    testLumpedInertiaQuad8FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaQuad8ReducedIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20ReducedIntegrationMatchesAnalyticValues,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
