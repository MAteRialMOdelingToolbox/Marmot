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
                                                               testPressureLoadTangent };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
