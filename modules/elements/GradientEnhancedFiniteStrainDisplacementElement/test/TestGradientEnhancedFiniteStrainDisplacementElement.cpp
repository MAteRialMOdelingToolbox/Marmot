#include "Marmot/GradientEnhancedFiniteStrainDisplacementElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 1.0 };

  const std::vector< double > quad8Coordinates =
    { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5 };

  const std::vector< double > hexa8Coordinates = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                                   0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };

  std::vector< double > hexa20Coordinates()
  {
    std::vector< double > c( hexa8Coordinates );
    const int             edges[12][2] = { { 0, 1 },
                                           { 1, 2 },
                                           { 2, 3 },
                                           { 3, 0 },
                                           { 4, 5 },
                                           { 5, 6 },
                                           { 6, 7 },
                                           { 7, 4 },
                                           { 0, 4 },
                                           { 1, 5 },
                                           { 2, 6 },
                                           { 3, 7 } };
    for ( const auto& e : edges )
      for ( int i = 0; i < 3; i++ )
        c.push_back( 0.5 * ( hexa8Coordinates[3 * e[0] + i] + hexa8Coordinates[3 * e[1] + i] ) );
    return c;
  }

  // An element made through the factory, with properties, state and initialization -- as a host does it.
  struct Setup {
    std::unique_ptr< MarmotElement > element;
    std::vector< double >            stateVars;
    std::vector< double >            elementProperties;

    Setup( const std::string& name, const std::vector< double >& coordinates )
      : element( MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) ), elementProperties( { 1.0 } )
    {
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty( ElementProperties( elementProperties.data(), elementProperties.size() ) );
      element->assignProperty(
        MarmotMaterialSection( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE", matProps.data(), matProps.size() ) );
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), stateVars.size() );
      element->initializeYourself();
      reset();
    }

    void reset()
    {
      std::fill( stateVars.begin(), stateVars.end(), 0.0 );
      element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    }

    int nDof() const { return element->getNDofPerElement(); }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& Q )
    {
      reset();
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof(), nDof() );
      element->computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0.0, 1.0 );
      return { P, K };
    }
  };

  // blocked element dof vector [u_1 .. u_n | N_1 .. N_n]: finite displacements, and a nonlocal field
  // above the damage threshold everywhere, so that the element is in the damaging (loading) branch
  Eigen::VectorXd damagingState( int nDim, int nNodes )
  {
    Eigen::VectorXd Q( nNodes * ( nDim + 1 ) );
    for ( int i = 0; i < nNodes * nDim; i++ )
      Q[i] = 0.04 * std::sin( 1.3 * ( i + 1 ) );
    for ( int A = 0; A < nNodes; A++ )
      Q[nNodes * nDim + A] = 3e-3 + 1e-3 * std::cos( 0.7 * ( A + 1 ) );
    return Q;
  }

  // central differences with a step per field: the nonlocal field lives on a much smaller scale (~1e-3)
  // than the displacements, and one common step would be dominated by truncation error there
  Eigen::MatrixXd numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                    const Eigen::VectorXd&                                            Q0,
                                    int                                                               nU )
  {
    Eigen::MatrixXd numK( f( Q0 ).size(), Q0.size() );
    for ( int j = 0; j < Q0.size(); j++ ) {
      const double    h  = j < nU ? 1e-7 : 1e-9;
      Eigen::VectorXd Qp = Q0, Qm = Q0;
      Qp[j] += h;
      Qm[j] -= h;
      numK.col( j ) = ( f( Qp ) - f( Qm ) ) / ( 2 * h );
    }
    return numK;
  }

  void checkNumericalTangent( const std::string& name, const std::vector< double >& coordinates, int nDim, int nNodes )
  {
    Setup s( name, coordinates );

    const Eigen::VectorXd Q0 = damagingState( nDim, nNodes );
    const auto [P, K]        = s.kernels( Q0 );

    const Eigen::MatrixXd numK = numericalTangent( [&]( const Eigen::VectorXd& Q ) { return s.kernels( Q ).first; },
                                                   Q0,
                                                   nNodes * nDim );

    const double relError = ( K - numK ).norm() / numK.norm();
    throwExceptionOnFailure( relError < 1e-7,
                             MakeString() << name << ": stiffness does not match the numerical tangent, relative error "
                                          << relError );

    // the coupling block dr_U/dN is small against the total, so check it on its own scale as well
    const int    nU       = nNodes * nDim;
    const double relErrUN = ( K.block( 0, nU, nU, nNodes ) - numK.block( 0, nU, nU, nNodes ) ).norm() /
                            numK.block( 0, nU, nU, nNodes ).norm();
    throwExceptionOnFailure( relErrUN < 1e-6,
                             MakeString()
                               << name << ": coupling block dr_U/dN inconsistent, relative error " << relErrUN );
  }

} // namespace

void testBasicProperties()
{
  Setup q8( "GCPE8UL", quad8Coordinates );
  throwExceptionOnFailure( q8.nDof() == 24, "GCPE8UL: 8 x (2 + 1) dofs expected" );
  throwExceptionOnFailure( q8.element->getNNodes() == 8 && q8.element->getNSpatialDimensions() == 2,
                           "GCPE8UL: 8 nodes in 2D" );
  throwExceptionOnFailure( q8.element->getElementShape() == "quad8", "GCPE8UL: shape" );
  throwExceptionOnFailure( q8.element->getNumberOfQuadraturePoints() == 9, "GCPE8UL: 3x3 Gauss points" );

  Setup q8r( "GCPE8RUL", quad8Coordinates );
  throwExceptionOnFailure( q8r.element->getNumberOfQuadraturePoints() == 4, "GCPE8RUL: 2x2 Gauss points" );

  Setup h8( "GC3D8UL", hexa8Coordinates );
  throwExceptionOnFailure( h8.nDof() == 32, "GC3D8UL: 8 x (3 + 1) dofs expected" );
  throwExceptionOnFailure( h8.element->getNumberOfQuadraturePoints() == 8, "GC3D8UL: 2x2x2 Gauss points" );
}

void testNodeFieldsAndPermutation()
{
  Setup      q8( "GCPE8UL", quad8Coordinates );
  const auto fields = q8.element->getNodeFields();
  throwExceptionOnFailure( fields.size() == 8, "one entry per node" );
  for ( const auto& f : fields )
    throwExceptionOnFailure( f == std::vector< std::string >{ "displacement", "nonlocal damage" }, "node fields" );

  // interleaved [u1x u1y N1 u2x ...] <- blocked [u | N]
  const auto pattern = q8.element->getDofIndicesPermutationPattern();
  for ( int A = 0; A < 8; A++ ) {
    throwExceptionOnFailure( pattern[2 * A] == 3 * A && pattern[2 * A + 1] == 3 * A + 1, "u permutation" );
    throwExceptionOnFailure( pattern[16 + A] == 3 * A + 2, "N permutation" );
  }
}

void testConsistentTangentPlaneStrainQuad8()
{
  checkNumericalTangent( "GCPE8UL", quad8Coordinates, 2, 8 );
}

void testConsistentTangentReducedQuad8()
{
  checkNumericalTangent( "GCPE8RUL", quad8Coordinates, 2, 8 );
}

void testConsistentTangentHexa8()
{
  checkNumericalTangent( "GC3D8UL", hexa8Coordinates, 3, 8 );
}

void testUndeformedStateIsStressFree()
{
  Setup      h8( "GC3D8UL", hexa8Coordinates );
  const auto P = h8.kernels( Eigen::VectorXd::Zero( h8.nDof() ) ).first;
  throwExceptionOnFailure( P.norm() < 1e-12, "the undeformed element must have zero residual" );
}

void testRigidRotationIsStressFree()
{
  // u = (R - I) X for a finite rotation about z, zero nonlocal field
  Setup           h8( "GC3D8UL", hexa8Coordinates );
  const double    phi = 0.6;
  Eigen::Matrix3d R;
  R << std::cos( phi ), -std::sin( phi ), 0, std::sin( phi ), std::cos( phi ), 0, 0, 0, 1;

  Eigen::VectorXd Q = Eigen::VectorXd::Zero( h8.nDof() );
  for ( int A = 0; A < 8; A++ ) {
    const Eigen::Vector3d X( &hexa8Coordinates[3 * A] );
    Q.segment< 3 >( 3 * A ) = ( R - Eigen::Matrix3d::Identity() ) * X;
  }

  const auto P = h8.kernels( Q ).first;
  throwExceptionOnFailure( P.norm() < 1e-9,
                           MakeString() << "rigid rotation must not load the element, |P| = " << P.norm() );
}

void testInternalForcesAreSelfEquilibrated()
{
  // translation invariance: the displacement residuals of every direction sum to zero
  Setup      h8( "GC3D8UL", hexa8Coordinates );
  const auto P = h8.kernels( damagingState( 3, 8 ) ).first;
  for ( int i = 0; i < 3; i++ ) {
    double sum = 0;
    for ( int A = 0; A < 8; A++ )
      sum += P[3 * A + i];
    throwExceptionOnFailure( std::abs( sum ) < 1e-10 * P.norm(), "displacement residuals must sum to zero" );
  }
}

void testPressureLoadTangent()
{
  Setup        q8( "GCPE8UL", quad8Coordinates );
  const double p    = 3.0;
  const int    face = 3;

  auto load = [&]( const Eigen::VectorXd& Q ) {
    Eigen::VectorXd P = Eigen::VectorXd::Zero( q8.nDof() );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( q8.nDof(), q8.nDof() );
    q8.element->computeDistributedLoad( MarmotElement::Pressure, P.data(), K.data(), face, &p, Q.data(), 0.0, 1.0 );
    return std::make_pair( P, K );
  };

  const Eigen::VectorXd Q0       = damagingState( 2, 8 );
  const auto [P0, K0]            = load( Q0 );
  const Eigen::MatrixXd numdP_dQ = numericalTangent( [&]( const Eigen::VectorXd& Q ) { return load( Q ).first; },
                                                     Q0,
                                                     16 );

  // hosts assemble R = P_int - P_ext, so a load reports K = -dP_ext/dQ (as DisplacementFiniteStrainULElement)
  throwExceptionOnFailure( P0.norm() > 0, "pressure must load the element" );
  throwExceptionOnFailure( ( K0 + numdP_dQ ).norm() < 1e-7 * numdP_dQ.norm(),
                           "pressure load tangent is inconsistent (expected K = -dP/dQ)" );
}

void testConsistentTangentHexa20()
{
  checkNumericalTangent( "GC3D20UL", hexa20Coordinates(), 3, 20 );
  checkNumericalTangent( "GC3D20RUL", hexa20Coordinates(), 3, 20 );
}

void testGeostaticInitialState()
{
  // constant hydrostatic in-situ stress: after the geostatic initialization, the undeformed element must carry it
  const double sigma           = -2.0;
  const double geostaticDef[6] = { sigma, 0.0, sigma, 1.0, 1.0, 1.0 }; // sigY(y=0), y1, sigY(y=1), y2, k11, k33

  for ( const auto& [name, coordinates, nDim, nNodes] : std::vector<
          std::tuple< std::string, std::vector< double >, int, int > >{ { "GCPE8UL", quad8Coordinates, 2, 8 },
                                                                        { "GC3D8UL", hexa8Coordinates, 3, 8 } } ) {
    Setup s( name, coordinates );
    s.element->setInitialConditions( MarmotElement::GeostaticStress, geostaticDef );

    Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
    Eigen::VectorXd Q = Eigen::VectorXd::Zero( s.nDof() );
    s.element->computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0.0, 1.0 );

    const double* tau = s.element->getStateView( "stress", 0 ).stateLocation;
    for ( int i = 0; i < 3; i++ )
      throwExceptionOnFailure( std::abs( tau[4 * i] - sigma ) < 1e-8 * std::abs( sigma ),
                               MakeString() << name << ": geostatic stress not reproduced, tau_" << i << i << " = "
                                            << tau[4 * i] );
  }
}

void testSurfaceTractionAndBodyForce()
{
  // resultants: traction x face area, body force x volume (unit square / cube, unit thickness)
  {
    Setup                 q8( "GCPE8UL", quad8Coordinates );
    const double          t[2] = { 0.7, -1.3 };
    Eigen::VectorXd       P    = Eigen::VectorXd::Zero( q8.nDof() );
    Eigen::MatrixXd       K    = Eigen::MatrixXd::Zero( q8.nDof(), q8.nDof() );
    const Eigen::VectorXd Q    = Eigen::VectorXd::Zero( q8.nDof() );
    q8.element->computeDistributedLoad( MarmotElement::SurfaceTraction, P.data(), K.data(), 1, t, Q.data(), 0., 1. );
    for ( int i = 0; i < 2; i++ ) {
      double sum = 0;
      for ( int A = 0; A < 8; A++ )
        sum += P[2 * A + i];
      throwExceptionOnFailure( std::abs( sum - t[i] ) < 1e-12, "surface traction resultant (Quad8)" );
    }

    const double b[2] = { 0.4, -2.1 };
    P.setZero();
    q8.element->computeBodyForce( P.data(), K.data(), b, Q.data(), 0., 1. );
    for ( int i = 0; i < 2; i++ ) {
      double sum = 0;
      for ( int A = 0; A < 8; A++ )
        sum += P[2 * A + i];
      throwExceptionOnFailure( std::abs( sum - b[i] ) < 1e-12, "body force resultant (Quad8)" );
    }
    throwExceptionOnFailure( P.tail( 8 ).norm() == 0.0, "loads act on the displacement field only" );
  }
  {
    Setup                 h8( "GC3D8UL", hexa8Coordinates );
    const double          b[3] = { 0.4, -2.1, 0.9 };
    Eigen::VectorXd       P    = Eigen::VectorXd::Zero( h8.nDof() );
    Eigen::MatrixXd       K    = Eigen::MatrixXd::Zero( h8.nDof(), h8.nDof() );
    const Eigen::VectorXd Q    = Eigen::VectorXd::Zero( h8.nDof() );
    h8.element->computeBodyForce( P.data(), K.data(), b, Q.data(), 0., 1. );
    for ( int i = 0; i < 3; i++ ) {
      double sum = 0;
      for ( int A = 0; A < 8; A++ )
        sum += P[3 * A + i];
      throwExceptionOnFailure( std::abs( sum - b[i] ) < 1e-12, "body force resultant (Hexa8)" );
    }
  }
}

void testPressureLoadTangentHexa8()
{
  Setup        h8( "GC3D8UL", hexa8Coordinates );
  const double p    = 3.0;
  const int    face = 2;

  auto load = [&]( const Eigen::VectorXd& Q ) {
    Eigen::VectorXd P = Eigen::VectorXd::Zero( h8.nDof() );
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( h8.nDof(), h8.nDof() );
    h8.element->computeDistributedLoad( MarmotElement::Pressure, P.data(), K.data(), face, &p, Q.data(), 0.0, 1.0 );
    return std::make_pair( P, K );
  };

  const Eigen::VectorXd Q0       = damagingState( 3, 8 );
  const auto [P0, K0]            = load( Q0 );
  const Eigen::MatrixXd numdP_dQ = numericalTangent( [&]( const Eigen::VectorXd& Q ) { return load( Q ).first; },
                                                     Q0,
                                                     24 );
  throwExceptionOnFailure( P0.norm() > 0, "pressure must load the element (Hexa8)" );
  throwExceptionOnFailure( ( K0 + numdP_dQ ).norm() < 1e-7 * numdP_dQ.norm(),
                           "pressure load tangent is inconsistent (Hexa8)" );
}

void testCoordinatesAndStateViews()
{
  Setup      q8( "GCPE8UL", quad8Coordinates );
  const auto center = q8.element->getCoordinatesAtCenter();
  throwExceptionOnFailure( std::abs( center[0] - 0.5 ) < 1e-14 && std::abs( center[1] - 0.5 ) < 1e-14,
                           "center of the unit square" );

  const auto qps = q8.element->getCoordinatesAtQuadraturePoints();
  throwExceptionOnFailure( (int)qps.size() == q8.element->getNumberOfQuadraturePoints(), "one coordinate per QP" );
  for ( const auto& x : qps )
    throwExceptionOnFailure( x[0] > 0 && x[0] < 1 && x[1] > 0 && x[1] < 1, "QPs inside the element" );

  // element state and material state are both reachable by name
  q8.kernels( damagingState( 2, 8 ) );
  throwExceptionOnFailure( q8.element->getStateView( "stress", 0 ).stateSize == 9, "element state: stress" );
  throwExceptionOnFailure( *q8.element->getStateView( "kappa", 0 ).stateLocation > 1e-3,
                           "material state: kappa beyond the damage threshold" );
}

void testUnsupportedRequestsThrow()
{
  Setup                 q8( "GCPE8UL", quad8Coordinates );
  const Eigen::VectorXd Q = Eigen::VectorXd::Zero( q8.nDof() );
  Eigen::VectorXd       P = Eigen::VectorXd::Zero( q8.nDof() );
  Eigen::MatrixXd       K = Eigen::MatrixXd::Zero( q8.nDof(), q8.nDof() );

  bool explicitThrew = false;
  try {
    q8.element->computeKernelsExplicit( Q.data(), Q.data(), P.data(), 0., 1. );
  }
  catch ( const std::runtime_error& ) {
    explicitThrew = true;
  }
  throwExceptionOnFailure( explicitThrew, "explicit kernels are not implemented and must say so" );

  bool loadThrew = false;
  try {
    const double load[1] = { 1.0 };
    q8.element->computeDistributedLoad( MarmotElement::SurfaceTorsion, P.data(), K.data(), 1, load, Q.data(), 0., 1. );
  }
  catch ( const std::invalid_argument& ) {
    loadThrew = true;
  }
  throwExceptionOnFailure( loadThrew, "an unsupported distributed load must be rejected" );

  bool initialConditionThrew = false;
  try {
    q8.element->setInitialConditions( MarmotElement::HydrostaticStress, nullptr );
  }
  catch ( const std::invalid_argument& ) {
    initialConditionThrew = true;
  }
  throwExceptionOnFailure( initialConditionThrew, "an unsupported initial condition must be rejected" );

  // plane stress is not registered, but the class accepts the section: it must refuse to compute
  using Element = GradientEnhancedFiniteStrainDisplacementElement< 2, 8 >;
  Element planeStress( 1, FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::PlaneStress );
  std::vector< double > elementProperties = { 1.0 };
  planeStress.assignNodeCoordinates( quad8Coordinates.data() );
  planeStress.assignProperty( ElementProperties( elementProperties.data(), elementProperties.size() ) );
  planeStress.assignProperty(
    MarmotMaterialSection( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE", matProps.data(), matProps.size() ) );
  std::vector< double > stateVars( planeStress.getNumberOfRequiredStateVars(), 0.0 );
  planeStress.assignStateVars( stateVars.data(), stateVars.size() );
  planeStress.initializeYourself();
  bool planeStressThrew = false;
  try {
    planeStress.computeKernels( Q.data(), Q.data(), P.data(), K.data(), 0., 1. );
  }
  catch ( const std::runtime_error& ) {
    planeStressThrew = true;
  }
  throwExceptionOnFailure( planeStressThrew, "plane stress must be refused" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testBasicProperties,
                                                               testNodeFieldsAndPermutation,
                                                               testConsistentTangentPlaneStrainQuad8,
                                                               testConsistentTangentReducedQuad8,
                                                               testConsistentTangentHexa8,
                                                               testUndeformedStateIsStressFree,
                                                               testRigidRotationIsStressFree,
                                                               testInternalForcesAreSelfEquilibrated,
                                                               testPressureLoadTangent,
                                                               testConsistentTangentHexa20,
                                                               testGeostaticInitialState,
                                                               testSurfaceTractionAndBodyForce,
                                                               testPressureLoadTangentHexa8,
                                                               testCoordinatesAndStateViews,
                                                               testUnsupportedRequestsThrow };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
