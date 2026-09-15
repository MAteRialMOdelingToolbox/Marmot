#include "Marmot/DisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::NumericalAlgorithms::Differentiation;

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

// ---------------------------------------------------------------------------------------------
// The setup below (nodal coordinates, material properties, state vars) is always kept alive in
// the SAME function scope as the element that uses it: assignNodeCoordinates()/assignProperty()
// store zero-copy Eigen::Map views (and MarmotMaterialFiniteStrain stores a raw pointer to its
// property array) directly into the caller-owned buffers rather than copying them, so factoring
// this setup into a helper that returns only the element -- while its backing vectors go out of
// scope -- leaves the element holding dangling pointers (confirmed via AddressSanitizer:
// heap-use-after-free when reading the node coordinates during computeKernels()).
// ---------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------
// Basic accessors
// ---------------------------------------------------------------------------------------------

void testBasicAccessorsHexa8()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  throwExceptionOnFailure( element->getNNodes() == 8, "Incorrect number of nodes." );
  throwExceptionOnFailure( element->getNSpatialDimensions() == 3, "Incorrect number of spatial dimensions." );
  throwExceptionOnFailure( element->getNDofPerElement() == 24, "Incorrect number of DOFs." );
  throwExceptionOnFailure( element->getElementShape() == "hexa8", "Incorrect element shape." );
  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == static_cast< int >( element->qps.size() ),
                           "getNumberOfQuadraturePoints() does not match the number of quadrature points." );

  const auto centerCoords = element->getCoordinatesAtCenter();
  throwExceptionOnFailure( centerCoords.size() == 3, "Incorrect dimension for center coordinates." );
  for ( int d = 0; d < 3; d++ )
    throwExceptionOnFailure( checkIfEqual( centerCoords[d], 0.5 ), "Incorrect center coordinate for unit cube." );

  const auto qpCoords = element->getCoordinatesAtQuadraturePoints();
  throwExceptionOnFailure( static_cast< int >( qpCoords.size() ) == element->getNumberOfQuadraturePoints(),
                           "getCoordinatesAtQuadraturePoints() returned the wrong number of entries." );
  for ( const auto& c : qpCoords )
    throwExceptionOnFailure( c.size() == 3, "Incorrect dimension for a quadrature point coordinate." );
}

// ---------------------------------------------------------------------------------------------
// computeKernels() / computeKernelsExplicit(): tangent consistency and implicit/explicit
// agreement.
//
// DisplacementFiniteStrainULElement has no independently derivable reference stiffness matrix
// (unlike the linear-elastic DisplacementFiniteElement), so its consistent tangent is instead
// verified the way AGENTS.md prescribes: against MarmotMathCore's numerical differentiation of
// the negative residual vector computeKernels() returns, at a moderately deformed configuration.
// ---------------------------------------------------------------------------------------------

void testComputeKernelsSolid3DTangentMatchesNumericalDifferentiation()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q0[i] = 0.02 * std::sin( 1.3 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd K_P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  K_P.setZero();
  K.setZero();
  element->computeKernels( Q0.data(), dQ.data(), K_P.data(), K.data(), 0.0, 1.0 );

  throwExceptionOnFailure( K.isApprox( K.transpose(), 1e-8 ),
                           "Hyperelastic stiffness matrix is not symmetric at a deformed configuration." );

  NumDiff::vector_to_vector_function_type residual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    Eigen::VectorXd Ptmp( nDof );
    Eigen::MatrixXd Ktmp( nDof, nDof );
    Ptmp.setZero();
    Ktmp.setZero();
    element->computeKernels( Q.data(), dQ.data(), Ptmp.data(), Ktmp.data(), 0.0, 1.0 );
    return Ptmp;
  };

  const Eigen::MatrixXd numK = NumDiff::centralDifference( residual, Q0 );

  throwExceptionOnFailure( checkIfEqual( K, numK, 1e-6 ),
                           "computeKernels() stiffness matrix does not match the numerical tangent (Solid, 3D)." );
}

void testComputeKernelsPlaneStrainTangentMatchesNumericalDifferentiation()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };         // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q0[i] = 0.02 * std::sin( 1.7 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  element->computeKernels( Q0.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );

  NumDiff::vector_to_vector_function_type residual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    Eigen::VectorXd Ptmp( nDof );
    Eigen::MatrixXd Ktmp( nDof, nDof );
    Ptmp.setZero();
    Ktmp.setZero();
    element->computeKernels( Q.data(), dQ.data(), Ptmp.data(), Ktmp.data(), 0.0, 1.0 );
    return Ptmp;
  };

  const Eigen::MatrixXd numK = NumDiff::centralDifference( residual, Q0 );

  throwExceptionOnFailure( checkIfEqual( K, numK, 1e-6 ),
                           "computeKernels() stiffness matrix does not match the numerical tangent "
                           "(PlaneStrain)." );
}

void testComputeKernelsPlaneStressThrowsNotImplemented()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStress;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 };
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int                   nDof = element->getNDofPerElement();
  const std::vector< double > Q( nDof, 0.0 );
  const std::vector< double > dQ( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );

  bool threw = false;
  try {
    element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeKernels() must throw for plane stress (not yet implemented)." );
}

void testComputeKernelsExplicitMatchesImplicitResidualSolid3D()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVecImplicit = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                        0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };
  const std::vector< double > nodeCoordsVecExplicit = nodeCoordsVecImplicit;
  const std::vector< double > matProps              = { 1.0, 1.0, 1.0 }; // K, G, density

  auto elementImplicit = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  elementImplicit->assignNodeCoordinates( nodeCoordsVecImplicit.data() );
  MarmotMaterialSection materialSectionImplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  elementImplicit->assignProperty( materialSectionImplicit );
  const int             nStateVarsTotalImplicit = elementImplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsImplicit( nStateVarsTotalImplicit, 0.0 );
  elementImplicit->assignStateVars( stateVarsImplicit.data(), nStateVarsTotalImplicit );
  elementImplicit->initializeYourself();

  auto elementExplicit = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  elementExplicit->assignNodeCoordinates( nodeCoordsVecExplicit.data() );
  MarmotMaterialSection materialSectionExplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  elementExplicit->assignProperty( materialSectionExplicit );
  const int             nStateVarsTotalExplicit = elementExplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsExplicit( nStateVarsTotalExplicit, 0.0 );
  elementExplicit->assignStateVars( stateVarsExplicit.data(), nStateVarsTotalExplicit );
  elementExplicit->initializeYourself();

  const int nDof = elementImplicit->getNDofPerElement();

  Eigen::VectorXd Q( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q[i] = 0.02 * std::sin( 0.9 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd PImplicit( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  PImplicit.setZero();
  K.setZero();
  elementImplicit->computeKernels( Q.data(), dQ.data(), PImplicit.data(), K.data(), 0.0, 1.0 );

  Eigen::VectorXd PExplicit( nDof );
  PExplicit.setZero();
  elementExplicit->computeKernelsExplicit( Q.data(), dQ.data(), PExplicit.data(), 0.0, 1.0 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( PImplicit ), Eigen::MatrixXd( PExplicit ), 1e-10 ),
                           "computeKernelsExplicit() residual does not match computeKernels() (Solid, 3D)." );
}

void testComputeKernelsExplicitMatchesImplicitResidualPlaneStrain()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  const std::vector< double > matProps      = { 1.0, 1.0, 1.0 }; // K, G, density
  const std::vector< double > elPropsVec    = { 1.0 };           // thickness

  auto elementImplicit = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  elementImplicit->assignNodeCoordinates( nodeCoordsVec.data() );
  MarmotMaterialSection materialSectionImplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  ElementProperties     elPropsImplicit( elPropsVec.data(), elPropsVec.size() );
  elementImplicit->assignProperty( elPropsImplicit );
  elementImplicit->assignProperty( materialSectionImplicit );
  const int             nStateVarsTotalImplicit = elementImplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsImplicit( nStateVarsTotalImplicit, 0.0 );
  elementImplicit->assignStateVars( stateVarsImplicit.data(), nStateVarsTotalImplicit );
  elementImplicit->initializeYourself();

  auto elementExplicit = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  elementExplicit->assignNodeCoordinates( nodeCoordsVec.data() );
  MarmotMaterialSection materialSectionExplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  ElementProperties     elPropsExplicit( elPropsVec.data(), elPropsVec.size() );
  elementExplicit->assignProperty( elPropsExplicit );
  elementExplicit->assignProperty( materialSectionExplicit );
  const int             nStateVarsTotalExplicit = elementExplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsExplicit( nStateVarsTotalExplicit, 0.0 );
  elementExplicit->assignStateVars( stateVarsExplicit.data(), nStateVarsTotalExplicit );
  elementExplicit->initializeYourself();

  const int nDof = elementImplicit->getNDofPerElement();

  Eigen::VectorXd Q( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q[i] = 0.02 * std::sin( 1.1 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd PImplicit( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  PImplicit.setZero();
  K.setZero();
  elementImplicit->computeKernels( Q.data(), dQ.data(), PImplicit.data(), K.data(), 0.0, 1.0 );

  Eigen::VectorXd PExplicit( nDof );
  PExplicit.setZero();
  elementExplicit->computeKernelsExplicit( Q.data(), dQ.data(), PExplicit.data(), 0.0, 1.0 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( PImplicit ), Eigen::MatrixXd( PExplicit ), 1e-10 ),
                           "computeKernelsExplicit() residual does not match computeKernels() (PlaneStrain)." );
}

// ---------------------------------------------------------------------------------------------
// setInitialConditions()
// ---------------------------------------------------------------------------------------------

void testSetInitialConditionsMaterialInitializationSetsUnitEigenDeformation()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  for ( const auto& qp : element->qps )
    throwExceptionOnFailure( checkIfEqual( qp.managedStateVars->F0_XX, 0.0 ), "F0_XX must start at zero." );

  element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

  for ( const auto& qp : element->qps ) {
    throwExceptionOnFailure( checkIfEqual( qp.managedStateVars->F0_XX, 1.0 ),
                             "setInitialConditions(MarmotMaterialInitialization) must set F0_XX to 1." );
    throwExceptionOnFailure( checkIfEqual( qp.managedStateVars->F0_YY, 1.0 ),
                             "setInitialConditions(MarmotMaterialInitialization) must set F0_YY to 1." );
    throwExceptionOnFailure( checkIfEqual( qp.managedStateVars->F0_ZZ, 1.0 ),
                             "setInitialConditions(MarmotMaterialInitialization) must set F0_ZZ to 1." );
  }
}

void testSetInitialConditionsGeostaticStressMatchesPrescribedStress()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };         // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  // The Newton iteration inside findEigenDeformationForEigenStress() starts from the qp's
  // *current* F0 (see setInitialConditions(GeostaticStress) below) and needs a non-degenerate
  // starting point; MarmotMaterialInitialization sets F0 to the identity for exactly this reason
  // and must run first, mirroring how a real analysis initializes an element.
  element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

  // Depth-independent (no vertical gradient) geostatic stress sigma_22 = -1 at y=0 and y=1, with
  // K0 = 1.0 in both horizontal directions: an isotropic (hydrostatic-like) target state sigma_11
  // = sigma_22 = sigma_33 = -1. Kept modest relative to K = G = 1: findEigenDeformationForEigenStress()
  // only budgets 5 Newton iterations, so a much larger target risks not converging in time.
  const std::vector< double > geostaticDefinition = { -1.0, 0.0, -1.0, 1.0, 1.0, 1.0 };

  element->setInitialConditions( MarmotElement::GeostaticStress, geostaticDefinition.data() );

  throwExceptionOnFailure( element->hasEigenDeformation,
                           "setInitialConditions(GeostaticStress) must set hasEigenDeformation to true." );

  const auto& qp0 = element->qps[0];

  Fastor::Tensor< double, 3, 3 > F0( 0.0 );
  F0( 0, 0 ) = qp0.managedStateVars->F0_XX;
  F0( 1, 1 ) = qp0.managedStateVars->F0_YY;
  F0( 2, 2 ) = qp0.managedStateVars->F0_ZZ;

  MarmotMaterialFiniteStrain::Deformation< 3 >          deformation{ F0 };
  MarmotMaterialFiniteStrain::TimeIncrement             timeIncrement{ 0.0, 1.0 };
  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response;
  response.stateVars = qp0.managedStateVars->materialStateVars.data();
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 > tangents;

  qp0.material->computeStress( response, tangents, deformation, timeIncrement );

  Fastor::Tensor< double, 3, 3 > expected( 0.0 );
  expected( 0, 0 ) = -1.0;
  expected( 1, 1 ) = -1.0;
  expected( 2, 2 ) = -1.0;

  throwExceptionOnFailure( checkIfEqual( response.tau, expected, 1e-3 ),
                           "The eigen deformation found by setInitialConditions(GeostaticStress) does not "
                           "reproduce the prescribed geostatic stress state." );
}

void testComputeKernelsPlaneStrainWithEigenDeformationTangentMatchesNumericalDifferentiation()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };         // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  // See the comment in testSetInitialConditionsGeostaticStressMatchesPrescribedStress(): the
  // eigen-deformation Newton iteration needs a non-degenerate starting F0.
  element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

  const std::vector< double > geostaticDefinition = { -1.0, 0.0, -1.0, 1.0, 0.5, 0.5 };
  element->setInitialConditions( MarmotElement::GeostaticStress, geostaticDefinition.data() );
  throwExceptionOnFailure( element->hasEigenDeformation, "hasEigenDeformation must be set before this test proceeds." );

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q0[i] = 0.01 * std::sin( 2.1 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  element->computeKernels( Q0.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );

  NumDiff::vector_to_vector_function_type residual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    Eigen::VectorXd Ptmp( nDof );
    Eigen::MatrixXd Ktmp( nDof, nDof );
    Ptmp.setZero();
    Ktmp.setZero();
    element->computeKernels( Q.data(), dQ.data(), Ptmp.data(), Ktmp.data(), 0.0, 1.0 );
    return Ptmp;
  };

  const Eigen::MatrixXd numK = NumDiff::centralDifference( residual, Q0 );

  throwExceptionOnFailure( checkIfEqual( K, numK, 1e-6 ),
                           "computeKernels() stiffness matrix does not match the numerical tangent "
                           "(PlaneStrain with eigen deformation)." );
}

// ---------------------------------------------------------------------------------------------
// computeConsistentInertia(), computeBodyForce(): both integrate a shape-function-weighted field
// over the element, so a uniform ("rigid-body") field in one direction must integrate to exactly
// that field times the corresponding total (mass or force), independent of element distortion,
// since the shape functions form a partition of unity.
// ---------------------------------------------------------------------------------------------

void testComputeConsistentInertiaConservesTotalMassPerDirection()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDim   = 3;
  const int nNodes = 8;
  const int nDof   = element->getNDofPerElement();

  double totalMassFromQuadrature = 0.0;
  for ( const auto& qp : element->qps )
    totalMassFromQuadrature += qp.material->getDensity( qp.managedStateVars->materialStateVars.data() ) * qp.J0xW;

  std::vector< double > M( nDof * nDof, 0.0 );
  element->computeConsistentInertia( M.data() );

  Eigen::Map< Eigen::MatrixXd > Mmat( M.data(), nDof, nDof );
  throwExceptionOnFailure( Mmat.isApprox( Mmat.transpose(), 1e-12 ), "Consistent mass matrix is not symmetric." );

  for ( int d = 0; d < nDim; d++ ) {
    Eigen::VectorXd uRigid = Eigen::VectorXd::Zero( nDof );
    for ( int a = 0; a < nNodes; a++ )
      uRigid[a * nDim + d] = 1.0;

    const double massInDirection = uRigid.transpose() * Mmat * uRigid;
    throwExceptionOnFailure( checkIfEqual( massInDirection, totalMassFromQuadrature, 1e-10 ),
                             "computeConsistentInertia() does not conserve the total element mass." );
  }
}

void testComputeBodyForceConservesTotalForcePerDirection()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDim   = 3;
  const int nNodes = 8;
  const int nDof   = element->getNDofPerElement();

  double totalVolume = 0.0;
  for ( const auto& qp : element->qps )
    totalVolume += qp.J0xW;

  const std::vector< double > load( { 1.0, 2.0, 3.0 } );
  const std::vector< double > QTotal( nDof, 0.0 );

  std::vector< double > P( nDof, 0.0 );
  std::vector< double > K( nDof * nDof, 0.0 ); // body force has zero stiffness contribution; unused here
  element->computeBodyForce( P.data(), K.data(), load.data(), QTotal.data(), 0.0, 1.0 );

  Eigen::Map< Eigen::VectorXd > Pvec( P.data(), nDof );

  for ( int d = 0; d < nDim; d++ ) {
    Eigen::VectorXd uRigid = Eigen::VectorXd::Zero( nDof );
    for ( int a = 0; a < nNodes; a++ )
      uRigid[a * nDim + d] = 1.0;

    const double forceInDirection = uRigid.dot( Pvec );
    throwExceptionOnFailure( checkIfEqual( forceInDirection, load[d] * totalVolume, 1e-10 ),
                             "computeBodyForce() does not integrate to the analytically expected total force." );
  }
}

// ---------------------------------------------------------------------------------------------
// computeCriticalTimeStepForExplicitDynamics()
// ---------------------------------------------------------------------------------------------

void testComputeCriticalTimeStepMatchesMaterialWaveSpeedForRegularHexa8()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int                   nDof = element->getNDofPerElement();
  const std::vector< double > QTotal( nDof, 0.0 );

  double criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  const auto& qp0 = element->qps[0];

  Fastor::Tensor< double, 3, 3 > F( 0.0 );
  F( 0, 0 ) = 1.0;
  F( 1, 1 ) = 1.0;
  F( 2, 2 ) = 1.0;

  const double c = qp0.material->getMaximumWaveSpeed( qp0.managedStateVars->materialStateVars.data(), F );

  // Hexa8 is a regular, linear (non-serendipity) element, so the mass-distribution factor (see
  // MarmotMassLumping.h) is exactly 1 for every node, and the characteristic length of the unit
  // cube is 2 * 0.5 = 1 (twice the smallest Jacobian singular value).
  const double expected = 1.0 / c;

  throwExceptionOnFailure( checkIfEqual( criticalTimeStep, expected, 1e-8 ),
                           "Hexa8 critical time step does not match the material wave speed for a unit mass "
                           "distribution factor." );
}

// ---------------------------------------------------------------------------------------------
// computeInternalEnergy(), getStateView(), computeDistributedLoad()
// ---------------------------------------------------------------------------------------------

void testComputeInternalEnergyMatchesSumOverQuadraturePoints()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q[i] = 0.02 * std::sin( 1.5 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );

  double expected = 0.0;
  for ( const auto& qp : element->qps )
    expected += qp.managedStateVars->totalStrainEnergy;

  double internalEnergy = 0.0;
  element->computeInternalEnergy( internalEnergy );

  throwExceptionOnFailure( checkIfEqual( internalEnergy, expected, 1e-12 ),
                           "computeInternalEnergy() does not match the sum over quadrature points." );
  throwExceptionOnFailure( internalEnergy > 0.0,
                           "computeInternalEnergy() should be strictly positive for a deformed hyperelastic "
                           "element." );
}

void testGetStateViewReturnsQuadraturePointManagedState()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 3, 8 >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 3, 8 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q[i] = 0.02 * std::sin( 0.5 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );

  const auto stressView = element->getStateView( "stress", 0 );
  throwExceptionOnFailure( stressView.stateSize == 9, "getStateView(\"stress\") returned the wrong size." );
  for ( int i = 0; i < 9; i++ )
    throwExceptionOnFailure( checkIfEqual( stressView.stateLocation[i], element->qps[0].managedStateVars->stress( i ) ),
                             "getStateView(\"stress\") does not point at the quadrature point's managed stress." );

  bool threwForUnknownName = false;
  try {
    element->getStateView( "this state does not exist", 0 );
  }
  catch ( const std::exception& ) {
    threwForUnknownName = true;
  }
  throwExceptionOnFailure( threwForUnknownName,
                           "getStateView() for an unrecognized state name must fall through to the material "
                           "and fail." );
}

void testComputeDistributedLoadThrowsForUnhandledLoadType()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteStrainULElement< 2, 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };         // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  const std::vector< double > load( 2, 0.0 );
  const std::vector< double > QTotal( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );

  bool threw = false;
  try {
    element->computeDistributedLoad( MarmotElement::SurfaceTorsion,
                                     P.data(),
                                     K.data(),
                                     0,
                                     load.data(),
                                     QTotal.data(),
                                     0.0,
                                     1.0 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeDistributedLoad() must throw for an unhandled load type." );
}

// ---------------------------------------------------------------------------------------------
// AxiSymmetricDisplacementFiniteStrainULElement: computeKernels()/computeKernelsExplicit() are
// overridden as *private* members there (matching MarmotElement's public virtual interface), so
// they are exercised through a MarmotElement* -- exactly how the element factory hands elements
// to callers, and the only way calling them compiles from outside the class.
// ---------------------------------------------------------------------------------------------

void testAxiSymmetricComputeKernelsTangentMatchesNumericalDifferentiation()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  // r in [1, 2], z in [0, 1]: kept away from the symmetry axis r = 0, since the axisymmetric
  // kernels divide by r.
  const std::vector< double > nodeCoordsVec = { 1.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 1.0 };

  std::unique_ptr< MarmotElement >
    element = std::make_unique< AxiSymmetricDisplacementFiniteStrainULElement< 4 > >( 1, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 1.0, 1.0, 1.0 }; // K, G, density
  MarmotMaterialSection       materialSection( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  element->assignProperty( materialSection );
  // AxiSymmetricDisplacementFiniteStrainULElement inherits initializeYourself() from the nDim=2
  // base, which unconditionally reads elementProperties[0] as a thickness (even though the
  // axisymmetric kernels themselves integrate 2*pi*r and never use it): a dummy value is
  // required here purely to satisfy that read.
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const int nDof = element->getNDofPerElement();

  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q0[i] = 0.01 * std::sin( 1.9 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  element->computeKernels( Q0.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );

  NumDiff::vector_to_vector_function_type residual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    Eigen::VectorXd Ptmp( nDof );
    Eigen::MatrixXd Ktmp( nDof, nDof );
    Ptmp.setZero();
    Ktmp.setZero();
    element->computeKernels( Q.data(), dQ.data(), Ptmp.data(), Ktmp.data(), 0.0, 1.0 );
    return Ptmp;
  };

  const Eigen::MatrixXd numK = NumDiff::centralDifference( residual, Q0 );

  throwExceptionOnFailure( checkIfEqual( K, numK, 1e-5 ),
                           "AxiSymmetricDisplacementFiniteStrainULElement::computeKernels() stiffness matrix "
                           "does not match the numerical tangent." );
}

void testAxiSymmetricComputeKernelsExplicitMatchesImplicitResidual()
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = DisplacementFiniteStrainULElement< 2, 4 >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 1.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 1.0 };
  const std::vector< double > matProps      = { 1.0, 1.0, 1.0 }; // K, G, density

  std::unique_ptr< MarmotElement >
    elementImplicit = std::make_unique< AxiSymmetricDisplacementFiniteStrainULElement< 4 > >( 1, intType, secType );
  elementImplicit->assignNodeCoordinates( nodeCoordsVec.data() );
  MarmotMaterialSection materialSectionImplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  elementImplicit->assignProperty( materialSectionImplicit );
  const std::vector< double > elPropsVecImplicit = { 1.0 }; // dummy thickness; see comment above
  ElementProperties           elPropsImplicit( elPropsVecImplicit.data(), elPropsVecImplicit.size() );
  elementImplicit->assignProperty( elPropsImplicit );
  const int             nStateVarsTotalImplicit = elementImplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsImplicit( nStateVarsTotalImplicit, 0.0 );
  elementImplicit->assignStateVars( stateVarsImplicit.data(), nStateVarsTotalImplicit );
  elementImplicit->initializeYourself();

  std::unique_ptr< MarmotElement >
    elementExplicit = std::make_unique< AxiSymmetricDisplacementFiniteStrainULElement< 4 > >( 1, intType, secType );
  elementExplicit->assignNodeCoordinates( nodeCoordsVec.data() );
  MarmotMaterialSection materialSectionExplicit( "COMPRESSIBLENEOHOOKE", matProps.data(), matProps.size() );
  elementExplicit->assignProperty( materialSectionExplicit );
  const std::vector< double > elPropsVecExplicit = { 1.0 }; // dummy thickness; see comment above
  ElementProperties           elPropsExplicit( elPropsVecExplicit.data(), elPropsVecExplicit.size() );
  elementExplicit->assignProperty( elPropsExplicit );
  const int             nStateVarsTotalExplicit = elementExplicit->getNumberOfRequiredStateVars();
  std::vector< double > stateVarsExplicit( nStateVarsTotalExplicit, 0.0 );
  elementExplicit->assignStateVars( stateVarsExplicit.data(), nStateVarsTotalExplicit );
  elementExplicit->initializeYourself();

  const int nDof = elementImplicit->getNDofPerElement();

  Eigen::VectorXd Q( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q[i] = 0.01 * std::sin( 0.8 * ( i + 1 ) );
  const Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof );

  Eigen::VectorXd PImplicit( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  PImplicit.setZero();
  K.setZero();
  elementImplicit->computeKernels( Q.data(), dQ.data(), PImplicit.data(), K.data(), 0.0, 1.0 );

  Eigen::VectorXd PExplicit( nDof );
  PExplicit.setZero();
  elementExplicit->computeKernelsExplicit( Q.data(), dQ.data(), PExplicit.data(), 0.0, 1.0 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( PImplicit ), Eigen::MatrixXd( PExplicit ), 1e-10 ),
                           "AxiSymmetricDisplacementFiniteStrainULElement::computeKernelsExplicit() residual "
                           "does not match computeKernels()." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testLumpedInertiaHexa8RegularElementIsPositiveAndConservesMass,
    testLumpedInertiaQuad8FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaQuad8ReducedIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20ReducedIntegrationMatchesAnalyticValues,
    testBasicAccessorsHexa8,
    testComputeKernelsSolid3DTangentMatchesNumericalDifferentiation,
    testComputeKernelsPlaneStrainTangentMatchesNumericalDifferentiation,
    testComputeKernelsPlaneStressThrowsNotImplemented,
    testComputeKernelsExplicitMatchesImplicitResidualSolid3D,
    testComputeKernelsExplicitMatchesImplicitResidualPlaneStrain,
    testSetInitialConditionsMaterialInitializationSetsUnitEigenDeformation,
    testSetInitialConditionsGeostaticStressMatchesPrescribedStress,
    testComputeKernelsPlaneStrainWithEigenDeformationTangentMatchesNumericalDifferentiation,
    testComputeConsistentInertiaConservesTotalMassPerDirection,
    testComputeBodyForceConservesTotalForcePerDirection,
    testComputeCriticalTimeStepMatchesMaterialWaveSpeedForRegularHexa8,
    testComputeInternalEnergyMatchesSumOverQuadraturePoints,
    testGetStateViewReturnsQuadraturePointManagedState,
    testComputeDistributedLoadThrowsForUnhandledLoadType,
    testAxiSymmetricComputeKernelsTangentMatchesNumericalDifferentiation,
    testAxiSymmetricComputeKernelsExplicitMatchesImplicitResidual,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
