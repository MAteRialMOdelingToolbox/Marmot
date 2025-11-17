#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

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
  const std::string                          matName  = "LINEARELASTIC"; // Linear Elastic
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
  const std::string                          matName  = "LINEARELASTIC"; // Linear Elastic
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
  double pNewDT      = 1.0; // Placeholder for adaptive time stepping (not used here)

  // Compute the stiffness matrix K and internal force vector P
  element->computeYourself( u.data(), dQ.data(), P.data(), K.data(), &currentTime, dt, pNewDT );

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
  const std::string                          matName  = "LINEARELASTIC"; // Linear Elastic
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

int main()
{
  auto tests = std::vector< std::function< void() > >{ testInstantiationAndBasicProperties,
                                                       testStiffnessMatrixCalculationPlaneStress,
                                                       testInitializeYourselfAndShapeFunctions };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
