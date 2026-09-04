#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <string>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

void testDefaultNamedPropertyInterface()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStress;

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );

  // callers (e.g. the EdelweissFE Cython wrapper) reach the named-property interface
  // through the MarmotElement base pointer, so exercise it the same way here
  MarmotElement* base = element.get();

  // an element that does not override the named-property interface must expose
  // no named properties, and assigning an unrecognized one must throw
  throwExceptionOnFailure( base->getPropertyNames().empty(),
                           "Default getPropertyNames() must be empty for an element without named properties." );

  const double dummyValue = 1.0;
  bool         threw      = false;
  try {
    base->assignProperty( "nonexistent property", &dummyValue );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "Default assignProperty(name, ...) must throw for an unrecognized named property." );
}

void testInstantiationAndBasicProperties()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  const auto                  secType       = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStress;

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const static std::vector< double > matProps = { 10000.0, 0.2, 1 };
  const std::string                  matName  = "LINEARELASTIC"; // Linear Elastic
  MarmotMaterialSection              materialSection( matName, matProps.data(), matProps.size() );
  const static std::vector< double > elPropsVec = { 0.1 };
  ElementProperties                  elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );

  element->initializeYourself();

  // Check Number of Quadrature Points (Gauss2x2 for Quad4 -> 4 QPs)
  int nQP = element->getNumberOfQuadraturePoints();
  throwExceptionOnFailure( nQP == 4, "Incorrect number of quadrature points." );

  // Check Element Shape
  std::string shape = element->getElementShape();
  throwExceptionOnFailure( shape == "quad4", "Incorrect element shape." );

  // Check Degrees of Freedom (nNodes * nDim = 4 * 2 = 8)
  int nDof = element->getNDofPerElement();
  throwExceptionOnFailure( nDof == 8, "Incorrect number of DOFs." );

  // Check Coordinates at Center (Should be 0.5, 0.5 for unit square)
  std::vector< double > centerCoords = element->getCoordinatesAtCenter();

  throwExceptionOnFailure( centerCoords.size() == nDim, "Incorrect dimension for center coordinates." );
  throwExceptionOnFailure( checkIfEqual( centerCoords[0], 0.5 ), "Incorrect X coordinate at center." );
  throwExceptionOnFailure( checkIfEqual( centerCoords[1], 0.5 ), "Incorrect Y coordinate at center." );
}

void testStiffnessMatrixCalculationPlaneStress()
{
  // Test setup:
  // Element Type: 2D Quad4
  // Integration: Full (2x2 Gauss)
  // Section Type: Plane Stress
  // Node Coordinates: (0,0), (6,0), (8,6), (2,6) - forms a general quadrilateral
  // Material: Linear Elastic (E=10000, nu=0.2, density=1)
  // Thickness: 1.0

  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;

  // Node coordinates defining a general quadrilateral
  const std::vector< double > nodeCoordsVec = { 0.0,
                                                0.0,   // Node 1 (x,y)
                                                6.0,
                                                0.0,   // Node 2 (x,y)
                                                8.0,
                                                6.0,   // Node 3 (x,y)
                                                2.0,
                                                6.0 }; // Node 4 (x,y)
  const auto                  secType       = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStress;

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  // Material properties: E = 10000, nu = 0.2, density = 1 (density not used in static analysis)
  const static std::vector< double > matProps = { 10000.0, 0.2, 1 };
  const std::string                  matName  = "LINEARELASTIC"; // Linear Elastic
  MarmotMaterialSection              materialSection( matName, matProps.data(), matProps.size() );

  // Element properties: thickness = 1.0
  const static std::vector< double > elPropsVec = { 1.0 };
  ElementProperties                  elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );

  element->initializeYourself();

  const int       nDof = element->getNDofPerElement();  // Should be nNodes * nDim = 4 * 2 = 8
  Eigen::VectorXd u    = Eigen::VectorXd::Zero( nDof ); // Nodal displacements (not used for K calculation directly)
  Eigen::VectorXd dQ   = Eigen::VectorXd::Zero( nDof ); // Incremental nodal displacements (not used for K calculation)
  Eigen::VectorXd P    = Eigen::VectorXd::Zero( nDof ); // Internal force vector
  Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( nDof, nDof ); // Stiffness matrix

  double currentTime = 0.0;
  double dt          = 1.0; // Dummy time step (not critical for linear elastic stiffness)
  // Compute the stiffness matrix K and internal force vector P
  element->computeKernels( u.data(), dQ.data(), P.data(), K.data(), currentTime, dt );

  // --- Stiffness Matrix Checks ---
  // The stiffness matrix K should be symmetric for linear elastic materials.
  double toleranceSymmetry = 1e-12; // Tolerance for symmetry check
  bool   isSymmetric       = K.isApprox( K.transpose(), toleranceSymmetry );
  throwExceptionOnFailure( isSymmetric, "Stiffness matrix is not symmetric." );

  // Reference stiffness matrix.
  // Values computed from an independent script for a plane-stress quad4 element
  // with nodes (0,0), (6,0), (8,6), (2,6), E=10000, nu=0.2, and thickness=1.0.
  Eigen::MatrixXd K_expected( nDof, nDof );
  // clang-format off
  K_expected << 4320.99,   868.06, -2932.1,    173.61, -1813.27, -1215.28,   424.38,   173.61,
                868.06,  3510.8,   1215.28,   -38.58, -1215.28,  -887.35,  -868.06, -2584.88,
               -2932.1,   1215.28,  5709.88, -2256.94,   424.38,  -868.06, -3202.16,  1909.72,
                173.61,   -38.58, -2256.94,  6983.02,   173.61, -2584.88,  1909.72, -4359.57,
               -1813.27, -1215.28,   424.38,   173.61,  4320.99,   868.06, -2932.1,    173.61,
               -1215.28,  -887.35,  -868.06, -2584.88,   868.06,  3510.8,   1215.28,   -38.58,
                424.38,  -868.06, -3202.16,  1909.72, -2932.1,   1215.28,  5709.88, -2256.94,
                173.61, -2584.88,  1909.72, -4359.57,   173.61,   -38.58, -2256.94,  6983.02;
  // clang-format on

  // Define a common tolerance for value checks.
  double valueTolerance = 1e-1;

  // Check the entire stiffness matrix against the expected values.
  throwExceptionOnFailure( checkIfEqual( K, K_expected, valueTolerance ),
                           "Stiffness matrix does not match expected values." );
}

void testInitializeYourselfAndShapeFunctions()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;

  // Unit square coordinates
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  const auto                  secType       = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStress;

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  // Dummy material properties (not strictly needed for this test but required by assignProperty)
  const static std::vector< double > matProps = { 10000.0, 0.2, 1 };
  const std::string                  matName  = "LINEARELASTIC"; // Linear Elastic
  MarmotMaterialSection              materialSection( matName, matProps.data(), matProps.size() );

  // Element properties: thickness = 1.0
  const static std::vector< double > elPropsVec = { 1.0 };
  ElementProperties                  elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection ); // Material assignment needed to avoid nullptrs if density is accessed

  // Assign dummy state variables
  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );

  element->initializeYourself();

  // --- Check properties of the first Quadrature Point ---
  // For a unit square element [0,1]x[0,1] mapped from canonical [-1,1]x[-1,1]:
  // Jacobian J = [0.5, 0; 0, 0.5], so detJ = 0.25.
  // For 2x2 Gauss quadrature, all weights are 1.0.
  // Thickness is 1.0. So, J0xW = detJ * weight * thickness = 0.25 * 1.0 * 1.0 = 0.25.

  const auto& qp0 = element->qps[0]; // Assuming consistent QP ordering

  throwExceptionOnFailure( checkIfEqual( qp0.detJ, 0.25 ), "Incorrect detJ for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.J0xW, 0.25 ), "Incorrect J0xW for QP0." );

  // --- Check B-matrix components for the first QP ---
  // Natural coordinates for the first Gauss point (see MarmotFiniteElement.h)
  const double xi_val  = 1.0 / std::sqrt( 3.0 );
  const double eta_val = 1.0 / std::sqrt( 3.0 );

  // Expected derivatives dN/dx, dN/dy for J_inv = [2,0; 0,2]
  // dNi/dx = 2 * dNi/dxi
  // dNi/dy = 2 * dNi/deta

  // For N1: dN1/dxi = -0.25*(1-eta), dN1/deta = -0.25*(1-xi)
  const double expected_dN1dx = 2.0 * ( -0.25 * ( 1.0 - eta_val ) ); // -0.5 * (1-eta_val)
  const double expected_dN1dy = 2.0 * ( -0.25 * ( 1.0 - xi_val ) );  // -0.5 * (1-xi_val)

  // B(0,0) = dN1/dx
  // B(1,1) = dN1/dy
  // B(2,0) = dN1/dy
  // B(2,1) = dN1/dx
  throwExceptionOnFailure( checkIfEqual( qp0.B( 0, 0 ), expected_dN1dx ), "Incorrect B(0,0) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 1, 1 ), expected_dN1dy ), "Incorrect B(1,1) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 2, 0 ), expected_dN1dy ), "Incorrect B(2,0) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 2, 1 ), expected_dN1dx ), "Incorrect B(2,1) for QP0." );

  // For N2: dN2/dxi = 0.25*(1-eta), dN2/deta = -0.25*(1+xi)
  const double expected_dN2dx = 2.0 * ( 0.25 * ( 1.0 - eta_val ) ); //  0.5 * (1-eta_val)
  const double expected_dN2dy = 2.0 * ( -0.25 * ( 1.0 + xi_val ) ); // -0.5 * (1+xi_val)

  // B(0,2) = dN2/dx
  // B(1,3) = dN2/dy
  // B(2,2) = dN2/dy
  // B(2,3) = dN2/dx
  throwExceptionOnFailure( checkIfEqual( qp0.B( 0, 2 ), expected_dN2dx ), "Incorrect B(0,2) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 1, 3 ), expected_dN2dy ), "Incorrect B(1,3) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 2, 2 ), expected_dN2dy ), "Incorrect B(2,2) for QP0." );
  throwExceptionOnFailure( checkIfEqual( qp0.B( 2, 3 ), expected_dN2dx ), "Incorrect B(2,3) for QP0." );
}

// ---------------------------------------------------------------------------------------------
// Lumped (diagonal) mass matrix tests
//
// computeLumpedInertia() uses the manifold-based scheme of Yang et al. (2017), mixing the
// high-order shape function N with the corresponding corner-node linear shape function N_lin
// via N_weighted = w*N + (1-w)*N_lin (only the corner entries receive the N_lin correction),
// with w = 1/2 by default and w = 1/3 special-cased for Hexa20 (see below). The tests below check
// the two properties any lumping scheme must satisfy to be usable in explicit dynamics: (1) no
// singular (zero or negative) nodal mass, and (2) conservation of the total element mass. Linear
// elements are included as a sanity baseline; quadratic elements (Quad8, Hexa20) get closer
// scrutiny, since the correction term is exactly what can drive a nodal mass towards (or
// through) zero.
// ---------------------------------------------------------------------------------------------

void testLumpedInertiaQuad4RegularElementIsPositiveAndConservesMass()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4 (linear)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 }; // unit square

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 10000.0, 0.2, density };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
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

  // Unit square, uniform density: by symmetry each corner carries exactly area*density/4.
  for ( int i = 0; i < nNodes; i++ ) {
    throwExceptionOnFailure( M[i * nDim] > 0.0, "Quad4 lumped mass entry is not strictly positive." );
    throwExceptionOnFailure( checkIfEqual( M[i * nDim], 0.25, 1e-12 ),
                             "Quad4 lumped mass entry does not match the analytic value." );
  }
}

void testLumpedInertiaHexa8RegularElementIsPositiveAndConservesMass()
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 8; // Hexa8 (linear)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::Solid;

  // Unit cube; node ordering per MarmotFiniteElement3D.cpp Hexa8::N.
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 };

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 10000.0, 0.2, density };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );

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
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStrain;

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

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 10000.0, 0.2, density };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  // Reference values from an independent symbolic (SymPy) double integration of
  // 0.5*N_serendipity + 0.5*N_bilinear over the reference square: corners = 1/12,
  // midsides = 1/6. Because the geometry map is affine here (straight edges, midsides at
  // exact midpoints), both full (3x3) and reduced (2x2) Gauss integrate the (low-degree)
  // integrand exactly, so both integration types must reproduce these same values.
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

void checkQuad8DistortedElementStaysNonSingular( FiniteElement::Quadrature::IntegrationTypes intType,
                                                 const std::string&                          label )
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 8;
  const int     elId    = 1;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStrain;

  // Curved-edge quarter-annulus sector (r_in=1, r_out=2, 0-90deg): a realistic curved-boundary
  // quadratic element, where the isoparametric mapping is no longer affine (unlike the regular
  // element above), exercising the mass-lumping scheme's robustness under distortion.
  const std::vector< double > nodeCoordsVec = {
    1.0,
    0.0,              // node 0: (r_in, 0deg)
    2.0,
    0.0,              // node 1: (r_out, 0deg)
    0.0,
    2.0,              // node 2: (r_out, 90deg)
    0.0,
    1.0,              // node 3: (r_in, 90deg)
    1.5,
    0.0,              // node 4: mid of edge 0-1
    std::sqrt( 2.0 ),
    std::sqrt( 2.0 ), // node 5: true arc midpoint of edge 1-2 (curved!)
    0.0,
    1.5,              // node 6: mid of edge 2-3
    std::sqrt( 0.5 ),
    std::sqrt( 0.5 )  // node 7: true arc midpoint of edge 3-0 (curved!)
  };

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 10000.0, 0.2, density };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  // Independently compute the total element mass from the same quadrature rule the element
  // itself uses (sum of rho * J0xW over all quadrature points), instead of relying on a
  // hardcoded reference value for this distorted (non-affine) geometry.
  double totalMassFromQuadrature = 0.0;
  for ( const auto& qp : element->qps ) {
    const double rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
    totalMassFromQuadrature += rho * qp.J0xW;
  }

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  double totalLumpedMass = 0.0;
  for ( int i = 0; i < nNodes; i++ ) {
    throwExceptionOnFailure( M[i * nDim] > 1e-8,
                             label + ": distorted Quad8 lumped mass entry is not strictly positive (node " +
                               std::to_string( i ) + ")." );
    totalLumpedMass += M[i * nDim];
  }

  throwExceptionOnFailure( checkIfEqual( totalLumpedMass, totalMassFromQuadrature, 1e-10 ),
                           label + ": distorted Quad8 lumped mass does not conserve the total element mass." );
}

void testLumpedInertiaQuad8DistortedElementFullIntegrationStaysNonSingular()
{
  checkQuad8DistortedElementStaysNonSingular( FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                              "FullIntegration" );
}

void testLumpedInertiaQuad8DistortedElementReducedIntegrationStaysNonSingular()
{
  checkQuad8DistortedElementStaysNonSingular( FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                                              "ReducedIntegration" );
}

void checkHexa20AnalyticLumpedMasses( FiniteElement::Quadrature::IntegrationTypes intType, const std::string& label )
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 20; // Hexa20 (quadratic serendipity)
  const int     elId    = 1;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::Solid;

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

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const double                density  = 1.0;
  const std::vector< double > matProps = { 10000.0, 0.2, density };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );

  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  std::vector< double > M( element->getNDofPerElement(), 0.0 );
  element->computeLumpedInertia( M.data() );

  // Reference values from an independent symbolic (SymPy) double integration of
  // 1/3*N_serendipity + 2/3*N_trilinear over the reference cube: corners = 1/24,
  // edges = 1/18. Both full (3x3x3) and reduced (2x2x2) Gauss integrate this exactly for a
  // regular (affinely-mapped) element, so both integration types must reproduce these values.
  //
  // REGRESSION GUARD: the default 1/2-1/2 split (used for every other element) makes the
  // negative corner contribution of the Hexa20 serendipity shape function EXACTLY cancel the
  // positive corner contribution of the trilinear shape function for any regular element: all 8
  // corners end up with EXACTLY ZERO lumped mass. Hexa20 is special-cased in computeLumpedInertia
  // to a 1/3-2/3 split for exactly this reason; this test checks both that the corner masses are
  // strictly positive and that they match the analytic reference for that special-cased split.
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
// Critical time step tests
//
// computeCriticalTimeStepForExplicitDynamics() scales the plain l/c estimate by the same
// mass-distribution factor the lumped-mass fractions above imply (see MarmotMassLumping.h):
// with n nodes and smallest lumped-mass fraction f_min, factor = sqrt(n * f_min). That gives
//
//   Quad4 / Hexa8 (linear): all fractions equal 1/n     -> factor = sqrt(n * 1/n)  = 1
//   Quad8:  corner = 1/12, midside = 1/6, n = 8         -> factor = sqrt(8/12)     = sqrt(2/3)
//   Hexa20: corner = 1/24, edge = 1/18, n = 20          -> factor = sqrt(20/24)    = sqrt(5/6)
//
// using the same analytic corner/edge fractions the lumped-mass tests above already establish.
// All four elements here share the same unit-size regular geometry and LINEARELASTIC material
// (E = 10000, nu = 0.2, density = 1), so the l/c part of the estimate is identical for all of
// them and the expected critical time step is just that common l/c times the factor above.
// ---------------------------------------------------------------------------------------------

double linearElasticLcEstimateForUnitRegularElement()
{
  // Reproduces the closed-form 3D stiffness diagonal getMaximumWaveSpeed() evaluates at zero
  // strain: C11 = E(1-nu) / ((1+nu)(1-2nu)) dominates the shear terms for nu = 0.2.
  constexpr double E       = 10000.0;
  constexpr double nu      = 0.2;
  constexpr double density = 1.0;
  const double     C11     = E * ( 1.0 - nu ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
  const double     c       = std::sqrt( C11 / density );
  constexpr double l       = 1.0; // unit-size regular element: 2 * smallest Jacobian singular value (0.5)
  return l / c;
}

void testCriticalTimeStepQuad4RegularElementMatchesUnitMassDistributionFactor()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4; // Quad4 (linear)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStrain;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 }; // unit square

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 10000.0, 0.2, 1.0 };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 }; // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const std::vector< double > QTotal( element->getNDofPerElement(), 0.0 );
  double                      criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  throwExceptionOnFailure( checkIfEqual( criticalTimeStep, linearElasticLcEstimateForUnitRegularElement(), 1e-8 ),
                           "Quad4 critical time step does not match the unit mass-distribution factor." );
}

void testCriticalTimeStepHexa8RegularElementMatchesUnitMassDistributionFactor()
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 8; // Hexa8 (linear)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::Solid;

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                                0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0 }; // unit cube

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 10000.0, 0.2, 1.0 };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );

  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const std::vector< double > QTotal( element->getNDofPerElement(), 0.0 );
  double                      criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  throwExceptionOnFailure( checkIfEqual( criticalTimeStep, linearElasticLcEstimateForUnitRegularElement(), 1e-8 ),
                           "Hexa8 critical time step does not match the unit mass-distribution factor." );
}

void testCriticalTimeStepQuad8RegularElementMatchesAnalyticMassDistributionFactor()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 8; // Quad8 (quadratic serendipity)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::PlaneStrain;

  // Unit square with midside nodes at exact edge midpoints (straight edges); same geometry as
  // checkQuad8AnalyticLumpedMasses().
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

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 10000.0, 0.2, 1.0 };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
  const std::vector< double > elPropsVec = { 1.0 }; // thickness
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );

  element->assignProperty( elProps );
  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const std::vector< double > QTotal( element->getNDofPerElement(), 0.0 );
  double                      criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  const double expected = linearElasticLcEstimateForUnitRegularElement() * std::sqrt( 2.0 / 3.0 );
  throwExceptionOnFailure( checkIfEqual( criticalTimeStep, expected, 1e-8 ),
                           "Quad8 critical time step does not match the analytic mass-distribution factor." );
}

void testCriticalTimeStepHexa20RegularElementMatchesAnalyticMassDistributionFactor()
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 20; // Hexa20 (quadratic serendipity)
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = DisplacementFiniteElement< nDim, nNodes >::SectionType::Solid;

  // Unit cube with edge-midside nodes at exact midpoints (straight edges); same geometry as
  // checkHexa20AnalyticLumpedMasses(). Node ordering per MarmotFiniteElement3D.cpp Hexa20::N:
  // 0-7 corners, 8-19 edge midsides.
  const std::vector< double > nodeCoordsVec = {
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, // bottom corners
    0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, // top corners
    0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.5, 0.0, // bottom-face edge midsides
    0.5, 0.0, 1.0, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0, 0.0, 0.5, 1.0, // top-face edge midsides
    0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 1.0, 0.5  // vertical edge midsides
  };

  auto element = std::make_unique< DisplacementFiniteElement< nDim, nNodes > >( elId, intType, secType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > matProps = { 10000.0, 0.2, 1.0 };
  MarmotMaterialSection       materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );

  element->assignProperty( materialSection );

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );
  element->initializeYourself();

  const std::vector< double > QTotal( element->getNDofPerElement(), 0.0 );
  double                      criticalTimeStep = 0.0;
  element->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );

  const double expected = linearElasticLcEstimateForUnitRegularElement() * std::sqrt( 5.0 / 6.0 );
  throwExceptionOnFailure( checkIfEqual( criticalTimeStep, expected, 1e-8 ),
                           "Hexa20 critical time step does not match the analytic mass-distribution factor." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testDefaultNamedPropertyInterface,
    testInstantiationAndBasicProperties,
    testStiffnessMatrixCalculationPlaneStress,
    testInitializeYourselfAndShapeFunctions,
    testLumpedInertiaQuad4RegularElementIsPositiveAndConservesMass,
    testLumpedInertiaHexa8RegularElementIsPositiveAndConservesMass,
    testLumpedInertiaQuad8FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaQuad8ReducedIntegrationMatchesAnalyticValues,
    testLumpedInertiaQuad8DistortedElementFullIntegrationStaysNonSingular,
    testLumpedInertiaQuad8DistortedElementReducedIntegrationStaysNonSingular,
    testLumpedInertiaHexa20FullIntegrationMatchesAnalyticValues,
    testLumpedInertiaHexa20ReducedIntegrationMatchesAnalyticValues,
    testCriticalTimeStepQuad4RegularElementMatchesUnitMassDistributionFactor,
    testCriticalTimeStepHexa8RegularElementMatchesUnitMassDistributionFactor,
    testCriticalTimeStepQuad8RegularElementMatchesAnalyticMassDistributionFactor,
    testCriticalTimeStepHexa20RegularElementMatchesAnalyticMassDistributionFactor,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
