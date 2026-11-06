#include "Marmot/EmbeddedBondElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Discretizations::ElementProperties;

namespace {

  /// A bond element created through the factory, together with the storage it needs
  struct Bond {
    std::unique_ptr< Marmot::Discretizations::Elements::Element > element;
    std::vector< double >                                         coordinates;
    std::vector< double >                                         properties;
    std::vector< double >                                         materialProperties;
    std::vector< double >                                         stateVars;
    int                                                           nDof, nDim;

    Bond( const std::string&    type,
          std::vector< double > coordinates_,
          std::vector< double > properties_,
          const std::string&    law,
          std::vector< double > lawProperties )
      : element( Marmot::Registry::ElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( std::move( properties_ ) ),
        materialProperties( std::move( lawProperties ) )
    {
      nDof = element->getNDofPerElement();
      nDim = element->getNSpatialDimensions();
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty( ElementProperties( properties.data(), static_cast< int >( properties.size() ) ) );
      element->assignProperty(
        Marmot::Discretizations::MaterialSection( law,
                                                  materialProperties.data(),
                                                  static_cast< int >( materialProperties.size() ) ) );
      element->initializeYourself();
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->setInitialConditions( Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization,
                                     nullptr );
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& U,
                                                           const Eigen::VectorXd& UOld,
                                                           bool                   commit )
    {
      const auto      backup = stateVars;
      Eigen::VectorXd P      = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K      = Eigen::MatrixXd::Zero( nDof, nDof );
      Eigen::VectorXd dU     = U - UOld;
      element->computeKernels( U.data(), dU.data(), P.data(), K.data(), 0.0, 1.0 );
      if ( !commit )
        stateVars = backup;
      return { P, K };
    }

    Eigen::MatrixXd numericalTangent( const Eigen::VectorXd& U, const Eigen::VectorXd& UOld, double h = 1e-8 )
    {
      Eigen::MatrixXd K( nDof, nDof );
      for ( int j = 0; j < nDof; j++ ) {
        Eigen::VectorXd Up = U, Um = U;
        Up( j ) += h;
        Um( j ) -= h;
        K.col( j ) = ( kernels( Up, UOld, false ).first - kernels( Um, UOld, false ).first ) / ( 2 * h );
      }
      return K;
    }

    /// nodal displacement vector of an affine displacement field u = a + H X
    Eigen::VectorXd affineField( const Eigen::VectorXd& a, const Eigen::MatrixXd& H ) const
    {
      Eigen::VectorXd U( nDof );
      for ( int i = 0; i < nDof / nDim; i++ )
        U.segment( i * nDim, nDim ) = a + H * Eigen::Map< const Eigen::VectorXd >( &coordinates[i * nDim], nDim );
      return U;
    }
  };

  // a distorted 8-node quadrilateral with curved edges (corners, then midsides) ...
  const std::vector< double > quad8 = { 0., 0., 4., -0.5, 4.5, 3., -0.5, 3.5, 2., -0.5, 4.5, 1.2, 2., 3.6, -0.4, 1.8 };
  // ... a distorted hexahedron
  const std::vector< double > hexa8 = { 0., 0., 0., 2., 0.,  0.2, 2.2, 2., 0.,  0.,  1.8, 0.,
                                        0., 0., 2., 2., 0.2, 2.,  2.,  2., 2.1, 0.1, 2.,  2. };

  std::vector< double > concat( std::vector< double > a, const std::vector< double >& b )
  {
    a.insert( a.end(), b.begin(), b.end() );
    return a;
  }

  const std::vector< double > linearBond = { 100., 1000. }; // Kt, Kn
  const double                tauMax     = 2.5 * std::sqrt( 30. );
  const std::vector< double > mc2010     = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., 1000. };

} // namespace

/// an affine displacement field shared by bar and host gives zero slip, for quadratic curved hosts (which requires the
/// inverse mapping to be correct) and in 3D
void testAffineFieldGivesNoSlip()
{
  {
    // inclined bar, ends inside the quad8
    Bond            bond( "EB2D2Q8",
               concat( { 0.5, 0.3, 3.5, 2.8 }, quad8 ),
                          { 1.0, -1.0, 1.0, 5 },
               "LINEARELASTICBONDSLIP",
               linearBond );
    Eigen::MatrixXd H( 2, 2 );
    H << 0.01, -0.03, 0.02, 0.005;
    auto [P, K] = bond.kernels( bond.affineField( Eigen::Vector2d( 0.3, -0.7 ), H ),
                                Eigen::VectorXd::Zero( bond.nDof ),
                                true );
    throwExceptionOnFailure( P.norm() < 1e-12,
                             Marmot::Utilities::MakeString()
                               << "EB2D2Q8: affine field must not slip, |P|=" << P.norm() );
  }
  {
    // a quadratic bar, curved, through the hexahedron
    Bond            bond( "EB3D3H8",
               concat( { 0.2, 0.3, 0.4, 1.7, 1.6, 1.5, 1.0, 0.8, 1.1 }, hexa8 ),
                          { 2.0, -1.0, 1.0 },
               "LINEARELASTICBONDSLIP",
               linearBond );
    Eigen::MatrixXd H = Eigen::MatrixXd::Random( 3, 3 ) * 0.01;
    auto [P, K]       = bond.kernels( bond.affineField( Eigen::Vector3d( 0.3, -0.7, 0.1 ), H ),
                                Eigen::VectorXd::Zero( bond.nDof ),
                                true );
    throwExceptionOnFailure( P.norm() < 1e-12,
                             Marmot::Utilities::MakeString()
                               << "EB3D3H8: affine field must not slip, |P|=" << P.norm() );
  }
}

/// pull-out of a bar from a fixed host: the bond force is the bond stress x perimeter x embedded length, along the bar;
/// a transverse displacement is resisted by the normal stiffness
void testPullOutForce()
{
  // bar along x from 0.5 to 3.5 at y = 1.5 inside the quad8, only the part [-0.5, 0.75] of the bar is in this host
  const double perimeter = 0.7;
  Bond         bond( "EB2D2Q8",
             concat( { 0.5, 1.5, 3.5, 1.5 }, quad8 ),
                     { perimeter, -0.5, 0.75, 3 },
             "LINEARELASTICBONDSLIP",
             linearBond );
  const double embeddedLength = 3.0 * ( 0.75 + 0.5 ) / 2;

  Eigen::VectorXd U = Eigen::VectorXd::Zero( bond.nDof );
  U( 0 ) = U( 2 ) = 0.01;  // axial slip
  U( 1 ) = U( 3 ) = 0.002; // transverse
  auto [P, K]     = bond.kernels( U, Eigen::VectorXd::Zero( bond.nDof ), true );

  const double axialForce = P( 0 ) + P( 2 );
  const double transverse = P( 1 ) + P( 3 );
  throwExceptionOnFailure( checkIfEqual( axialForce, 100. * 0.01 * perimeter * embeddedLength, 1e-12 ),
                           Marmot::Utilities::MakeString() << "pull-out force " << axialForce );
  throwExceptionOnFailure( checkIfEqual( transverse, 1000. * 0.002 * perimeter * embeddedLength, 1e-12 ),
                           Marmot::Utilities::MakeString() << "transverse force " << transverse );
  // equilibrium: the host takes the opposite force
  throwExceptionOnFailure( checkIfEqual( P.segment( 4, 16 ).reshaped( 2, 8 ).rowwise().sum().norm(),
                                         std::hypot( axialForce, transverse ),
                                         1e-12 ),
                           "host reaction" );
  throwExceptionOnFailure( P.sum() < 1e-12 && P.sum() > -1e-12, "equilibrium" );
  throwExceptionOnFailure( checkIfEqual( *bond.element->getStateView( "slip", 0 ).stateLocation, 0.01, 1e-14 ),
                           "slip state" );
}

/// consistent tangent with a path dependent, softening bond-slip law, for an inclined bar in a distorted hexahedron
void testTangentModelCode2010()
{
  Bond bond( "EB3D2H8",
             concat( { 0.3, 0.2, 0.4, 1.6, 1.7, 1.5 }, hexa8 ),
             { 1.5, -1.0, 1.0, 4 },
             "MODELCODE2010BONDSLIP",
             mc2010 );

  const Eigen::Vector3d t    = Eigen::Vector3d( 1.3, 1.5, 1.1 ).normalized();
  Eigen::VectorXd       UOld = Eigen::VectorXd::Zero( bond.nDof );
  // pull the bar out to the softening branch in steps
  for ( int i = 1; i <= 5; i++ ) {
    Eigen::VectorXd U   = Eigen::VectorXd::Zero( bond.nDof );
    U.segment< 3 >( 0 ) = t * 0.6 * i;
    U.segment< 3 >( 3 ) = t * 0.6 * i;
    bond.kernels( U, UOld, true );
    UOld = U;
  }
  Eigen::VectorXd U = UOld;
  U.segment< 3 >( 0 ) += t * 0.2 + Eigen::Vector3d( 0.01, -0.02, 0.0 );
  U.segment< 3 >( 3 ) += t * 0.1;
  U.segment( 6, 6 ) += Eigen::VectorXd::Constant( 6, 0.005 );

  const auto [P, K] = bond.kernels( U, UOld, false );
  const auto numK   = bond.numericalTangent( U, UOld );
  throwExceptionOnFailure( checkIfEqual< double >( K, numK, 1e-5 * numK.cwiseAbs().maxCoeff() ),
                           "bond tangent does not match the numerical one" );
}

void testBarOutsideOfHostThrows()
{
  bool threw = false;
  try {
    Bond bond( "EB2D2Q8",
               concat( { 0.5, 1.5, 9.5, 1.5 }, quad8 ),
               { 1.0, -1.0, 1.0 },
               "LINEARELASTICBONDSLIP",
               linearBond );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "a bond point outside of the host must be rejected" );
}

/// for tiny displacements, the large-slip element (channel in a quad8, window of one bar element) gives the forces of
/// the small-slip element

/// for tiny displacements, the large-slip element (channel in a quad8, window of one bar element) gives the forces of
/// the small-slip element
void testLargeSlipEqualsSmallSlipForSmallDisplacements()
{
  const std::vector< double > bar = { 0.5, 0.3, 3.5, 2.8 };
  Bond small( "EB2D2Q8", concat( bar, quad8 ), { 1.0, -1.0, 1.0, 3 }, "LINEARELASTICBONDSLIP", linearBond );
  Bond large( "EBLS2D2Q8W1",
              concat( quad8, bar ),
              { 1.0, 0.0, 1.0, 1.0, 1.0, 3 },
              "LINEARELASTICBONDSLIP",
              linearBond );

  Eigen::VectorXd USmall = Eigen::VectorXd::Random( small.nDof ) * 1e-6;
  // the same nodal displacements in the node order of the large-slip element: host, then bar
  Eigen::VectorXd ULarge( large.nDof );
  ULarge << USmall.segment( 4, 16 ), USmall.segment( 0, 4 );

  const auto [PSmall, KSmall] = small.kernels( USmall, Eigen::VectorXd::Zero( small.nDof ), false );
  const auto [PLarge, KLarge] = large.kernels( ULarge, Eigen::VectorXd::Zero( large.nDof ), false );
  Eigen::VectorXd PLargeReordered( small.nDof );
  PLargeReordered << PLarge.segment( 16, 4 ), PLarge.segment( 0, 16 );
  throwExceptionOnFailure( checkIfEqual< double >( PLargeReordered, PSmall, 1e-6 * PSmall.cwiseAbs().maxCoeff() ),
                           "large-slip forces must equal the small-slip ones for small displacements" );
  throwExceptionOnFailure( std::abs( KLarge.sum() - KSmall.sum() ) < 1e-4 * KSmall.cwiseAbs().sum(),
                           "large-slip tangent must equal the small-slip one for small displacements" );
}

/// a bar slid by 1.0 along its axis out of a fixed host: the channel point the bar end has passed loses its bond, the
/// others carry Kt x 1.0
void testLargeSlipPullOutLosesBond()
{
  const std::vector< double > bar = { 0.5, 1.5, 3.5, 1.5 };
  Bond                        large( "EBLS2D2Q8W1",
              concat( quad8, bar ),
                                     { 1.0, 0.0, 1.0, 1.0, 1.0, 3 },
              "LINEARELASTICBONDSLIP",
              linearBond );
  Eigen::VectorXd             U = Eigen::VectorXd::Zero( large.nDof );
  U( 16 ) = U( 18 ) = 1.0;
  const auto [P, K] = large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), true );

  // Lobatto points at x = 0.5, 2.0, 3.5 with tributary parts [0, 0.5], [0.5, 2.5], [2.5, 3] of the channel; the bar
  // start slid to 1.0: covered fractions 0, 0.75, 1, smoothed by 3x^2 - 2x^3 to 0, 0.84375, 1
  const double expected = 100. * 1.0 * ( 2.0 * 0.84375 + 0.5 * 1.0 );
  throwExceptionOnFailure( checkIfEqual( P( 16 ) + P( 18 ), expected, 1e-9 ),
                           Marmot::Utilities::MakeString() << "bar force " << P( 16 ) + P( 18 ) << " != " << expected );
  throwExceptionOnFailure( checkIfEqual( *large.element->getStateView( "active", 0 ).stateLocation, 0.0, 0.0 ) &&
                             checkIfEqual( *large.element->getStateView( "active", 2 ).stateLocation, 1.0, 0.0 ),
                           "the channel point the bar end has passed must be inactive" );
  throwExceptionOnFailure( checkIfEqual( *large.element->getStateView( "slip", 2 ).stateLocation, 1.0, 1e-12 ),
                           "slip of the active channel point" );
  throwExceptionOnFailure( checkIfEqual( *large.element->getStateView( "covered fraction", 1 ).stateLocation,
                                         0.75,
                                         1e-12 ),
                           "covered fraction of the middle channel point" );
}

void testLargeSlipWindowExceededThrows()
{
  const std::vector< double > bar = { 0.5, 1.5, 3.5, 1.5 };
  // the window start is not the end of the bar: sliding past it must be reported
  Bond large( "EBLS2D2Q8W1", concat( quad8, bar ), { 1.0, 0.0, 1.0, 0.0, 1.0 }, "LINEARELASTICBONDSLIP", linearBond );
  Eigen::VectorXd U = Eigen::VectorXd::Zero( large.nDof );
  U( 16 ) = U( 18 ) = 1.0;
  bool threw        = false;
  try {
    large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), false );
  }
  catch ( const std::runtime_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "a slip beyond the window must throw" );
}

/// the (internally finite-difference) tangent against an independent central difference, with a softening law and a
/// window of 3 quadratic bar elements in a hexahedron
void testLargeSlipTangent3D()
{
  // a straight quadratic bar of 3 elements through the hexahedron, chain order e0 m0 e1 m1 e2 m2 e3
  std::vector< double > chain;
  const Eigen::Vector3d a( 0.3, 0.2, 0.4 ), b( 1.6, 1.7, 1.5 );
  for ( int i = 0; i <= 6; i++ ) {
    const Eigen::Vector3d X = a + ( b - a ) * i / 6.0;
    chain.insert( chain.end(), { X( 0 ), X( 1 ), X( 2 ) } );
  }
  Bond large( "EBLS3D3H8W3", concat( hexa8, chain ), { 1.5, 0.4, 2.6, 1.0, 1.0, 4 }, "MODELCODE2010BONDSLIP", mc2010 );
  const Eigen::Vector3d t    = ( b - a ).normalized();
  Eigen::VectorXd       UOld = Eigen::VectorXd::Zero( large.nDof );
  for ( int i = 1; i <= 3; i++ ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( large.nDof );
    for ( int n = 0; n < 7; n++ )
      U.segment< 3 >( 24 + 3 * n ) = t * 0.1 * i;
    large.kernels( U, UOld, true );
    UOld = U;
  }
  Eigen::VectorXd U = UOld;
  for ( int n = 0; n < 7; n++ )
    U.segment< 3 >( 24 + 3 * n ) += t * 0.05 + Eigen::Vector3d( 0.002, -0.001, 0.0005 ) * n;
  U.segment( 0, 24 ) += Eigen::VectorXd::Constant( 24, 0.001 );

  const auto [P, K] = large.kernels( U, UOld, false );
  const auto numK   = large.numericalTangent( U, UOld, 1e-6 );
  throwExceptionOnFailure( checkIfEqual< double >( K, numK, 1e-4 * numK.cwiseAbs().maxCoeff() ),
                           "large-slip tangent does not match an independent numerical one" );
}

/// bar elements of different lengths: a channel point at the end of the channel (Gauss-Lobatto) belongs to its own
/// bar element, not to the (longer) next one; pushing the bar back by d gives Kt d p L of the channel exactly
void testLargeSlipUnequalBarElements()
{
  // chain of 2 bar elements: 0.5 -> 1.5 (length 1, the channel) -> 3.5 (length 2)
  const std::vector< double > chain = { 0.5, 1.5, 1.5, 1.5, 3.5, 1.5 };
  Bond large( "EBLS2D2Q8W2", concat( quad8, chain ), { 1.0, 0.0, 1.0, 1.0, 1.0 }, "LINEARELASTICBONDSLIP", linearBond );
  Eigen::VectorXd U = Eigen::VectorXd::Zero( large.nDof );
  const double    d = 1e-3;
  U( 16 ) = U( 18 ) = U( 20 ) = -d; // pushed back: no part of the channel is left by the bar
  const auto [P, K]           = large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), false );
  const double barForce       = P( 16 ) + P( 18 ) + P( 20 );
  throwExceptionOnFailure( checkIfEqual( barForce, -100. * d * 1.0 * 1.0, 1e-12 ),
                           Marmot::Utilities::MakeString()
                             << "bar force " << barForce << " != Kt d p L = " << -100. * d );
}

void testFactoryNames()
{
  for ( const auto& name : {
          "EB2D2Q4",
          "EB2D2Q8",
          "EB2D3Q4",
          "EB2D3Q8",
          "EB3D2T4",
          "EB3D2T10",
          "EB3D2H8",
          "EB3D2H20",
          "EB3D3T4",
          "EB3D3T10",
          "EB3D3H8",
          "EB3D3H20",
          "EBLS2D2Q4W1",
          "EBLS2D3Q8W48",
          "EBLS3D2H8W7",
          "EBLS3D3H20W12",
        } ) {
    auto el = std::unique_ptr< Marmot::Discretizations::Elements::Element >(
      Marmot::Registry::ElementFactory::createElement( name, 1 ) );
    throwExceptionOnFailure( el != nullptr, Marmot::Utilities::MakeString() << name << " must be registered" );
  }
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactoryNames,
    testAffineFieldGivesNoSlip,
    testPullOutForce,
    testTangentModelCode2010,
    testBarOutsideOfHostThrows,
    testLargeSlipEqualsSmallSlipForSmallDisplacements,
    testLargeSlipPullOutLosesBond,
    testLargeSlipWindowExceededThrows,
    testLargeSlipTangent3D,
    testLargeSlipUnequalBarElements,
  } );
  return 0;
}
