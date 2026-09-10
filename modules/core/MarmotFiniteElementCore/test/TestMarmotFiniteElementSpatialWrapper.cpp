#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::NumericalAlgorithms::Differentiation;

// ---------------------------------------------------------------------------------------------
// All tests below embed a 2-node 1D truss (DisplacementFiniteElement<1,2>, UniaxialStress) into
// 2D ambient space, exactly matching the production "T2D2" element registered in
// DisplacementFiniteElementRegistration.cpp: MarmotElementSpatialWrapper(nDim=2, nDimChild=1,
// nNodes=2, nRhsChild=2, rhsIndicesToBeWrapped={0,1}, nIndicesToBeWrapped=2, child).
// ---------------------------------------------------------------------------------------------

namespace {

  std::unique_ptr< MarmotElementSpatialWrapper > makeT2D2( const std::vector< double >& nodeCoordsVec,
                                                           const std::vector< double >& matProps,
                                                           std::vector< double >&       stateVarsOut )
  {
    auto child = std::make_unique<
      DisplacementFiniteElement< 1, 2 > >( 1,
                                           FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                           DisplacementFiniteElement< 1, 2 >::SectionType::UniaxialStress );

    static constexpr int indicesToBeWrapped[] = { 0, 1 };

    auto wrapper = std::make_unique< MarmotElementSpatialWrapper >( 2,
                                                                    1,
                                                                    2,
                                                                    2,
                                                                    indicesToBeWrapped,
                                                                    2,
                                                                    std::move( child ) );

    wrapper->assignNodeCoordinates( nodeCoordsVec.data() );

    MarmotMaterialSection materialSection( "LINEARELASTIC", matProps.data(), matProps.size() );
    // assignProperty(ElementProperties&) stores a zero-copy view into this array (see
    // DisplacementFiniteElement::assignProperty), so it must outlive the element; `static`
    // gives it program lifetime rather than dangling once this factory function returns.
    static const std::vector< double > elPropsVec = { 1.0 }; // cross-section area
    ElementProperties                  elProps( elPropsVec.data(), elPropsVec.size() );
    wrapper->assignProperty( elProps );
    wrapper->assignProperty( materialSection );

    const int nStateVarsTotal = wrapper->getNumberOfRequiredStateVars();
    stateVarsOut.assign( nStateVarsTotal, 0.0 );
    wrapper->assignStateVars( stateVarsOut.data(), nStateVarsTotal );
    wrapper->initializeYourself();

    return wrapper;
  }

} // namespace

// ---------------------------------------------------------------------------------------------
// Simple pass-through accessors
// ---------------------------------------------------------------------------------------------

void testBasicAccessorsDelegateToChild()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 }; // length-5 truss
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  throwExceptionOnFailure( wrapper->getNNodes() == 2, "getNNodes() must delegate to the wrapper's own nNodes." );
  throwExceptionOnFailure( wrapper->getNSpatialDimensions() == 2,
                           "getNSpatialDimensions() must return the ambient dimension." );
  throwExceptionOnFailure( wrapper->getNDofPerElement() == 4,
                           "getNDofPerElement() must return the ambient (unprojected) DOF count." );
  throwExceptionOnFailure( wrapper->getElementShape() == "bar2", "getElementShape() must delegate to the child." );
  throwExceptionOnFailure( wrapper->getNumberOfQuadraturePoints() > 0,
                           "getNumberOfQuadraturePoints() must delegate to the child." );

  const auto nodeFields = wrapper->getNodeFields();
  throwExceptionOnFailure( static_cast< int >( nodeFields.size() ) == 2,
                           "getNodeFields() must delegate to the child." );

  throwExceptionOnFailure( wrapper->getPropertyNames().empty(),
                           "getPropertyNames() must delegate to the child (empty by default)." );

  bool threw = false;
  try {
    const double dummy = 1.0;
    wrapper->assignProperty( "nonexistent property", &dummy );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "assignProperty(name, ...) must delegate to the child and throw." );
}

void testDofIndicesPermutationPatternExpandsChildPattern()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const auto pattern = wrapper->getDofIndicesPermutationPattern();
  throwExceptionOnFailure( static_cast< int >( pattern.size() ) == 4,
                           "getDofIndicesPermutationPattern() must have one entry per ambient DOF." );

  // Every projected (both, here) child DOF expands into nDim=2 consecutive ambient entries; the
  // pattern must at least be a valid permutation of 0..3 (no repeats, all indices in range).
  std::vector< int > sorted = pattern;
  std::sort( sorted.begin(), sorted.end() );
  for ( int i = 0; i < 4; i++ )
    throwExceptionOnFailure( sorted[i] == i, "getDofIndicesPermutationPattern() must be a permutation of 0..nDof-1." );
}

// ---------------------------------------------------------------------------------------------
// computeKernels(): tangent consistency, using the same "reset quadrature-point state, apply the
// full vector as dQ=Q from a zero baseline" convention as the incremental hypoelastic child
// element itself requires (see TestDisplacementFiniteElement.cpp).
// ---------------------------------------------------------------------------------------------

void testComputeKernelsTangentMatchesNumericalDifferentiation()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int  nDof       = wrapper->getNDofPerElement();
  const auto resetState = [&]() { std::fill( stateVars.begin(), stateVars.end(), 0.0 ); };

  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nDof; i++ )
    Q0[i] = 0.01 * std::sin( 1.3 * ( i + 1 ) );

  resetState();
  Eigen::VectorXd P( nDof );
  Eigen::MatrixXd K( nDof, nDof );
  P.setZero();
  K.setZero();
  wrapper->computeKernels( Q0.data(), Q0.data(), P.data(), K.data(), 0.0, 1.0 );

  NumDiff::vector_to_vector_function_type residual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    resetState();
    Eigen::VectorXd Ptmp( nDof );
    Eigen::MatrixXd Ktmp( nDof, nDof );
    Ptmp.setZero();
    Ktmp.setZero();
    wrapper->computeKernels( Q.data(), Q.data(), Ptmp.data(), Ktmp.data(), 0.0, 1.0 );
    return Ptmp;
  };

  const Eigen::MatrixXd numK = NumDiff::centralDifference( residual, Q0 );

  throwExceptionOnFailure( checkIfEqual( K, numK, 1e-6 ),
                           "MarmotElementSpatialWrapper::computeKernels() stiffness matrix does not match "
                           "the numerical tangent." );
}

void testComputeKernelsMatchesAnalyticAxialStiffnessForAxisAlignedTruss()
{
  // A truss along the x-axis reduces to a textbook 1D bar: K = EA/L * [[1,-1],[-1,1]] in the
  // AXIAL direction, and zero stiffness transverse to it (since UniaxialStress ignores lateral
  // displacement entirely).
  constexpr double            length        = 2.0;
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, length, 0.0 };
  constexpr double            E = 10000.0, A = 1.0;
  const std::vector< double > matProps = { E, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int                   nDof = wrapper->getNDofPerElement();
  const std::vector< double > Q( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );
  wrapper->computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0.0, 1.0 );

  Eigen::Map< Eigen::MatrixXd > Kmat( K.data(), nDof, nDof );

  const double axialStiffness = E * A / length;
  // DOF layout: [node0_x, node0_y, node1_x, node1_y]
  throwExceptionOnFailure( checkIfEqual( Kmat( 0, 0 ), axialStiffness, 1e-6 ),
                           "Axial stiffness K(0,0) does not match EA/L for an axis-aligned truss." );
  throwExceptionOnFailure( checkIfEqual( Kmat( 2, 2 ), axialStiffness, 1e-6 ),
                           "Axial stiffness K(2,2) does not match EA/L for an axis-aligned truss." );
  throwExceptionOnFailure( checkIfEqual( Kmat( 0, 2 ), -axialStiffness, 1e-6 ),
                           "Axial stiffness K(0,2) does not match -EA/L for an axis-aligned truss." );
  throwExceptionOnFailure( checkIfEqual( Kmat( 1, 1 ), 0.0, 1e-8 ),
                           "Transverse stiffness K(1,1) must be zero for a UniaxialStress truss." );
}

void testComputeKernelsExplicitThrowsNotImplemented()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int                   nDof = wrapper->getNDofPerElement();
  const std::vector< double > Q( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );

  bool threw = false;
  try {
    wrapper->computeKernelsExplicit( Q.data(), Q.data(), P.data(), 0.0, 1.0 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computeKernelsExplicit() is not overridden by the wrapper, so it must fall back to "
                           "MarmotElement's default (throwing) implementation." );
}

// ---------------------------------------------------------------------------------------------
// computeBodyForce(): a body force applied in the ambient frame, projected onto the truss axis
// and back, must still integrate to the ambient-frame total force times the truss length (the
// wrapper's P/P^T round trip should be exact for a linear, single-child-element case).
// ---------------------------------------------------------------------------------------------

void testComputeBodyForceConservesTotalForce()
{
  constexpr double            length        = 2.0;
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, length, 0.0 }; // axis-aligned truss
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int                   nDof = wrapper->getNDofPerElement();
  const std::vector< double > load = { 3.0, 0.0 }; // force per unit length, along the truss axis
  const std::vector< double > QTotal( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );

  wrapper->computeBodyForce( P.data(), K.data(), load.data(), QTotal.data(), 0.0, 1.0 );

  // Cross-section area = 1, so total axial force = load[0] * length.
  const double totalForceX = P[0] + P[2];
  const double totalForceY = P[1] + P[3];
  throwExceptionOnFailure( checkIfEqual( totalForceX, load[0] * length, 1e-8 ),
                           "computeBodyForce() through the wrapper does not integrate to load * length in x." );
  throwExceptionOnFailure( checkIfEqual( totalForceY, 0.0, 1e-8 ),
                           "computeBodyForce() through the wrapper introduced a spurious transverse force." );
}

// ---------------------------------------------------------------------------------------------
// setInitialConditions(): must delegate to the child.
// ---------------------------------------------------------------------------------------------

void testSetInitialConditionsDelegatesToChild()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  // DisplacementFiniteElement<1,...>::setInitialConditions only handles
  // MarmotMaterialInitialization and MarmotMaterialStateVars (throws) and rejects GeostaticStress
  // implicitly via its default branch, since nDim=1 there. MarmotMaterialInitialization must not
  // throw.
  wrapper->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

  bool threw = false;
  try {
    wrapper->setInitialConditions( MarmotElement::MarmotMaterialStateVars, nullptr );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "setInitialConditions() must delegate to the child and propagate its throw." );
}

// ---------------------------------------------------------------------------------------------
// getStateView(): "MarmotElementSpatialWrapper.T" is handled directly by the wrapper; anything
// else must delegate to the child.
// ---------------------------------------------------------------------------------------------

void testGetStateViewReturnsTransformationMatrixAndDelegatesOtherwise()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 }; // direction (0.6, 0.8)
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const auto tView = wrapper->getStateView( "MarmotElementSpatialWrapper.T", 0 );
  throwExceptionOnFailure( tView.stateSize == 2, "getStateView(\"...T\") returned the wrong size." );
  throwExceptionOnFailure( checkIfEqual( tView.stateLocation[0], 0.6, 1e-10 ) &&
                             checkIfEqual( tView.stateLocation[1], 0.8, 1e-10 ),
                           "getStateView(\"...T\") does not return the truss's direction cosines." );

  // Any other name must delegate to the child (e.g. "stress", which DisplacementFiniteElement
  // manages itself).
  const int                   nDof = wrapper->getNDofPerElement();
  const std::vector< double > Q( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );
  wrapper->computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0.0, 1.0 );

  const auto stressView = wrapper->getStateView( "stress", 0 );
  throwExceptionOnFailure( stressView.stateSize == 6, "getStateView(\"stress\") did not delegate to the child." );
}

// ---------------------------------------------------------------------------------------------
// getCoordinatesAtCenter() / getCoordinatesAtQuadraturePoints()
// ---------------------------------------------------------------------------------------------

void testGetCoordinatesAtCenterMatchesChildCoordinatesInAmbientSpace()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 }; // length-5 truss
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const auto center = wrapper->getCoordinatesAtCenter();

  // The truss midpoint in ambient space is the average of its two endpoints: (1.5, 2.0).
  throwExceptionOnFailure( static_cast< int >( center.size() ) == 2,
                           "getCoordinatesAtCenter() returned the wrong dimension." );
  throwExceptionOnFailure( checkIfEqual( center[0], 1.5, 1e-10 ) && checkIfEqual( center[1], 2.0, 1e-10 ),
                           "getCoordinatesAtCenter() does not return the truss midpoint in ambient space." );
}

void testGetCoordinatesAtQuadraturePointsMatchesChildCoordinatesInAmbientSpace()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const auto qpCoords = wrapper->getCoordinatesAtQuadraturePoints();
  throwExceptionOnFailure( static_cast< int >( qpCoords.size() ) == wrapper->getNumberOfQuadraturePoints(),
                           "getCoordinatesAtQuadraturePoints() returned the wrong number of entries." );

  for ( const auto& coords : qpCoords ) {
    throwExceptionOnFailure( static_cast< int >( coords.size() ) == 2,
                             "getCoordinatesAtQuadraturePoints() returned the wrong dimension." );
    // Every quadrature point must lie exactly on the line from (0,0) to (3,4), i.e. y/x = 4/3
    // (equivalently 4*x - 3*y == 0), and within the segment's bounding box.
    throwExceptionOnFailure( checkIfEqual( 4.0 * coords[0] - 3.0 * coords[1], 0.0, 1e-8 ),
                             "getCoordinatesAtQuadraturePoints() point does not lie on the truss axis." );
    throwExceptionOnFailure( coords[0] >= -1e-8 && coords[0] <= 3.0 + 1e-8,
                             "getCoordinatesAtQuadraturePoints() point lies outside the truss segment." );
  }
}

// ---------------------------------------------------------------------------------------------
// computeDistributedLoad(): forwards to the child (with load/QTotal in the correct order --
// verified by inspection against MarmotElement's declared signature, since the only production
// child of this wrapper is a 1D bar, and DisplacementFiniteElement<1,...>::computeDistributedLoad
// always throws before it would ever read those arguments -- see
// FiniteElement::BoundaryElement's constructor, which has no 1D "Bar2" parent-shape branch).
// ---------------------------------------------------------------------------------------------

void testComputeDistributedLoadDelegatesAndThrowsForUnsupported1DChild()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int                   nDof = wrapper->getNDofPerElement();
  const std::vector< double > load( 1, 1.0 );
  const std::vector< double > QTotal( nDof, 0.0 );
  std::vector< double >       P( nDof, 0.0 );
  std::vector< double >       K( nDof * nDof, 0.0 );

  bool threw = false;
  try {
    wrapper
      ->computeDistributedLoad( MarmotElement::Pressure, P.data(), K.data(), 0, load.data(), QTotal.data(), 0.0, 1.0 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computeDistributedLoad() must delegate to the child and propagate its throw (a 1D "
                           "child has no boundary-element support)." );
}

// ---------------------------------------------------------------------------------------------
// MarmotElementSpatialWrapper does not override computeLumpedInertia(), computeConsistentInertia(),
// computeCriticalTimeStepForExplicitDynamics(), or computeInternalEnergy(), so calling any of them
// falls through to MarmotElement's default (throwing) implementation -- otherwise untested,
// since every OTHER element in this codebase does override all four.
// ---------------------------------------------------------------------------------------------

void testUnoverriddenOptionalMethodsFallBackToBaseClassDefaults()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 3.0, 4.0 };
  const std::vector< double > matProps      = { 10000.0, 0.2, 1.0 };
  std::vector< double >       stateVars;

  auto wrapper = makeT2D2( nodeCoordsVec, matProps, stateVars );

  const int nDof = wrapper->getNDofPerElement();

  bool threw = false;
  try {
    std::vector< double > M( nDof, 0.0 );
    wrapper->computeLumpedInertia( M.data() );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeLumpedInertia() must fall back to the base class default (throwing)." );

  threw = false;
  try {
    std::vector< double > M( nDof * nDof, 0.0 );
    wrapper->computeConsistentInertia( M.data() );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeConsistentInertia() must fall back to the base class default (throwing)." );

  threw = false;
  try {
    double                      criticalTimeStep = 0.0;
    const std::vector< double > QTotal( nDof, 0.0 );
    wrapper->computeCriticalTimeStepForExplicitDynamics( criticalTimeStep, QTotal.data() );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "computeCriticalTimeStepForExplicitDynamics() must fall back to the base class "
                           "default (throwing)." );

  threw = false;
  try {
    double internalEnergy = 0.0;
    wrapper->computeInternalEnergy( internalEnergy );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "computeInternalEnergy() must fall back to the base class default (throwing)." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testBasicAccessorsDelegateToChild,
    testDofIndicesPermutationPatternExpandsChildPattern,
    testComputeKernelsTangentMatchesNumericalDifferentiation,
    testComputeKernelsMatchesAnalyticAxialStiffnessForAxisAlignedTruss,
    testComputeKernelsExplicitThrowsNotImplemented,
    testComputeBodyForceConservesTotalForce,
    testSetInitialConditionsDelegatesToChild,
    testGetStateViewReturnsTransformationMatrixAndDelegatesOtherwise,
    testGetCoordinatesAtCenterMatchesChildCoordinatesInAmbientSpace,
    testGetCoordinatesAtQuadraturePointsMatchesChildCoordinatesInAmbientSpace,
    testComputeDistributedLoadDelegatesAndThrowsForUnsupported1DChild,
    testUnoverriddenOptionalMethodsFallBackToBaseClassDefaults,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
