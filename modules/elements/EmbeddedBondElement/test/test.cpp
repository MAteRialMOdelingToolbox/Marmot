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
  // with Coulomb friction mu = 0.4 and crushing at 90 MPa
  const std::vector< double > mc2010Friction = { tauMax, 1.0, 2.0, 10.0, 0.4, 0.4 * tauMax, 200., 1000., 0.4, 90. };

  std::vector< double > shifted( std::vector< double > coordinates, int nDim, const std::vector< double >& offset )
  {
    for ( size_t i = 0; i < coordinates.size(); i++ )
      coordinates[i] += offset[i % nDim];
    return coordinates;
  }

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

/// far from the origin (offset 1e5 times the element size), the host parent coordinates, the forces and the tangent
/// are those of the element at the origin: the inverse mapping works on coordinates relative to the element
void testOffsetCoordinates()
{
  const auto coordinates3D = concat( { 0.3, 0.2, 0.4, 1.6, 1.7, 1.5 }, hexa8 );
  const auto coordinates2D = concat( { 0.5, 1.5, 3.5, 1.5 }, quad8 );
  for ( const double offset : { 1e5, -3e5 } ) {
    Bond            a( "EB3D2H8", coordinates3D, { 1.5, -1.0, 1.0, 4 }, "MODELCODE2010BONDSLIP", mc2010 );
    Bond            b( "EB3D2H8",
            shifted( coordinates3D, 3, { offset, -0.7 * offset, 0.3 * offset } ),
                       { 1.5, -1.0, 1.0, 4 },
            "MODELCODE2010BONDSLIP",
            mc2010 );
    Eigen::VectorXd U   = Eigen::VectorXd::Zero( a.nDof );
    U.segment< 3 >( 0 ) = Eigen::Vector3d( 0.3, 0.2, 0.1 );
    U.segment< 3 >( 3 ) = Eigen::Vector3d( 0.25, 0.3, 0.1 );
    const auto [Pa, Ka] = a.kernels( U, Eigen::VectorXd::Zero( a.nDof ), false );
    const auto [Pb, Kb] = b.kernels( U, Eigen::VectorXd::Zero( b.nDof ), false );
    throwExceptionOnFailure( checkIfEqual< double >( Pa, Pb, 1e-8 * Pa.norm() ) &&
                               checkIfEqual< double >( Ka, Kb, 1e-8 * Ka.norm() ),
                             Marmot::Utilities::MakeString() << "EB3D2H8 at offset " << offset << ": |dP| = "
                                                             << ( Pa - Pb ).norm() << ", |P| = " << Pa.norm() );

    Bond            c( "EB2D2Q8", coordinates2D, { 0.7, -0.5, 0.75, 3 }, "LINEARELASTICBONDSLIP", linearBond );
    Bond            d( "EB2D2Q8",
            shifted( coordinates2D, 2, { offset, 0.4 * offset } ),
                       { 0.7, -0.5, 0.75, 3 },
            "LINEARELASTICBONDSLIP",
            linearBond );
    Eigen::VectorXd V = Eigen::VectorXd::Zero( c.nDof );
    V( 0 ) = V( 2 ) = 0.01;
    V( 1 )          = 0.002;
    const auto Pc   = c.kernels( V, Eigen::VectorXd::Zero( c.nDof ), false ).first;
    const auto Pd   = d.kernels( V, Eigen::VectorXd::Zero( d.nDof ), false ).first;
    throwExceptionOnFailure( checkIfEqual< double >( Pc, Pd, 1e-8 * Pc.norm() ),
                             Marmot::Utilities::MakeString() << "EB2D2Q8 at offset " << offset );
  }
}

/// a bar in a linear triangle (CPS3/CPE3 host): an affine field gives no slip, a pull-out gives the closed-form force
void testTriangleHost()
{
  const std::vector< double > tri3 = { 0., 0., 4., 0.5, 1., 3. };
  Bond                        bond( "EB2D2T3",
             concat( { 0.8, 0.6, 2.8, 0.9 }, tri3 ),
                                    { 0.5, -1.0, 1.0 },
             "LINEARELASTICBONDSLIP",
             linearBond );
  Eigen::MatrixXd             H( 2, 2 );
  H << 0.01, -0.02, 0.03, 0.005;
  const auto P = bond
                   .kernels( bond.affineField( Eigen::Vector2d( 0.1, -0.2 ), H ),
                             Eigen::VectorXd::Zero( bond.nDof ),
                             false )
                   .first;
  throwExceptionOnFailure( P.norm() < 1e-12,
                           Marmot::Utilities::MakeString() << "EB2D2T3: affine field slips, |P| = " << P.norm() );

  const Eigen::Vector2d t        = Eigen::Vector2d( 2.0, 0.3 ).normalized();
  Eigen::VectorXd       U        = Eigen::VectorXd::Zero( bond.nDof );
  U.segment< 2 >( 0 )            = 0.01 * t;
  U.segment< 2 >( 2 )            = 0.01 * t;
  const auto            PPull    = bond.kernels( U, Eigen::VectorXd::Zero( bond.nDof ), false ).first;
  const double          length   = std::hypot( 2.0, 0.3 );
  const Eigen::Vector2d barForce = PPull.segment< 2 >( 0 ) + PPull.segment< 2 >( 2 );
  throwExceptionOnFailure( checkIfEqual< double >( barForce, Eigen::Vector2d( 100. * 0.01 * 0.5 * length * t ), 1e-12 ),
                           "EB2D2T3: pull-out force" );
}

/// invalid element properties are rejected, and so is a change of the number of integration points after the bond
/// law was assigned (it would destroy the laws)
void testInvalidElementProperties()
{
  const auto coordinates = concat( { 0.5, 1.5, 3.5, 1.5 }, quad8 );
  for ( const auto& properties : std::vector< std::vector< double > >{ { 0.0, -0.5, 0.5 },
                                                                       { 1.0, 0.5, -0.5 },
                                                                       { 1.0, -1.5, 0.5 },
                                                                       { 1.0, -0.5, 1.2 },
                                                                       { 1.0, -0.5, 0.5, 7 },
                                                                       { 1.0, -0.5 } } ) {
    bool threw = false;
    try {
      Bond bond( "EB2D2Q8", coordinates, properties, "LINEARELASTICBONDSLIP", linearBond );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw,
                             Marmot::Utilities::MakeString() << "element properties [" << properties[0] << ", "
                                                             << properties[1] << ", ...] must be rejected" );
  }

  auto element = std::unique_ptr< Marmot::Discretizations::Elements::Element >(
    Marmot::Registry::ElementFactory::createElement( "EB2D2Q8", 1 ) );
  element->assignNodeCoordinates( coordinates.data() );
  element->assignProperty( Marmot::Discretizations::MaterialSection( "LINEARELASTICBONDSLIP", linearBond.data(), 2 ) );
  const std::vector< double > properties = { 1.0, -0.5, 0.5, 4 };
  bool                        threw      = false;
  try {
    element->assignProperty( ElementProperties( properties.data(), 4 ) );
  }
  catch ( const std::logic_error& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "changing the integration points after the bond law must be rejected" );
}

/// consistent tangent with Coulomb friction, with and without the normal cap: the coupling of the tangential bond
/// stress to the normal slip enters the element tangent. The last increment pushes the bar further out normally, so
/// the points on the cap stay clearly on it (an increment tangential to the cap would put the trial pressure on the
/// kink pMax, where central differences average the elastic and the cap branch).
void testTangentWithFriction()
{
  const auto withoutCap = std::vector< double >( mc2010Friction.begin(), mc2010Friction.end() - 1 );
  for ( const auto& lawProperties : { mc2010Friction, withoutCap } ) {
    Bond bond( "EB3D2H8",
               concat( { 0.3, 0.2, 0.4, 1.6, 1.7, 1.5 }, hexa8 ),
               { 1.5, -1.0, 1.0, 3 },
               "MODELCODE2010BONDSLIP",
               lawProperties );

    const Eigen::Vector3d t    = Eigen::Vector3d( 1.3, 1.5, 1.1 ).normalized();
    const Eigen::Vector3d n    = t.cross( Eigen::Vector3d::UnitX() ).normalized();
    Eigen::VectorXd       UOld = Eigen::VectorXd::Zero( bond.nDof );
    for ( int i = 1; i <= 4; i++ ) {
      Eigen::VectorXd U   = Eigen::VectorXd::Zero( bond.nDof );
      U.segment< 3 >( 0 ) = t * 0.6 * i + n * 0.02 * i;
      U.segment< 3 >( 3 ) = t * 0.6 * i + n * 0.03 * i;
      bond.kernels( U, UOld, true );
      UOld = U;
    }
    Eigen::VectorXd U = UOld;
    U.segment< 3 >( 0 ) += t * 0.2 + n * 0.02;
    U.segment< 3 >( 3 ) += t * 0.1 + n * 0.03;

    const auto [P, K] = bond.kernels( U, UOld, false );
    const auto numK   = bond.numericalTangent( U, UOld );
    throwExceptionOnFailure( checkIfEqual< double >( K, numK, 1e-5 * numK.cwiseAbs().maxCoeff() ),
                             Marmot::Utilities::MakeString() << "bond tangent with friction (" << lawProperties.size()
                                                             << " properties) does not match the numerical one" );
    // without the cap, the pressure grows with the normal slip and the friction makes the tangent unsymmetric; on the
    // cap the pressure is constant and there is no coupling
    if ( lawProperties.size() == withoutCap.size() )
      throwExceptionOnFailure( !K.isApprox( K.transpose() ), "the friction coupling makes the tangent unsymmetric" );
  }
}

/// for small displacements, the large-slip element (channel in a quad8, window of one bar element) gives the forces and
/// the tangent of the small-slip element, entry by entry (the node orders differ: host first here)
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
  Eigen::PermutationMatrix< Eigen::Dynamic > largeToSmall( small.nDof );
  for ( int j = 0; j < 16; j++ )
    largeToSmall.indices()( j ) = 4 + j;
  for ( int j = 0; j < 4; j++ )
    largeToSmall.indices()( 16 + j ) = j;
  const Eigen::VectorXd ULarge = largeToSmall.transpose() * USmall;

  const auto [PSmall, KSmall]              = small.kernels( USmall, Eigen::VectorXd::Zero( small.nDof ), false );
  const auto [PLarge, KLarge]              = large.kernels( ULarge, Eigen::VectorXd::Zero( large.nDof ), false );
  const Eigen::VectorXd PLargeInSmallOrder = largeToSmall * PLarge;
  const Eigen::MatrixXd KLargeInSmallOrder = largeToSmall * KLarge * largeToSmall.transpose();
  throwExceptionOnFailure( checkIfEqual< double >( PLargeInSmallOrder, PSmall, 1e-5 * PSmall.cwiseAbs().maxCoeff() ),
                           "large-slip forces must equal the small-slip ones for small displacements" );
  throwExceptionOnFailure( checkIfEqual< double >( KLargeInSmallOrder, KSmall, 1e-5 * KSmall.cwiseAbs().maxCoeff() ),
                           Marmot::Utilities::MakeString()
                             << "large-slip tangent must equal the small-slip one entry by entry, max difference "
                             << ( KLargeInSmallOrder - KSmall ).cwiseAbs().maxCoeff() );
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

/// a straight quadratic bar of 3 elements through the hexahedron, window order e0 m0 e1 m1 e2 m2 e3
std::vector< double > straightWindow3D( const Eigen::Vector3d& a, const Eigen::Vector3d& b, int nNodes )
{
  std::vector< double > window;
  for ( int i = 0; i < nNodes; i++ ) {
    const Eigen::Vector3d X = a + ( b - a ) * i / double( nNodes - 1 );
    window.insert( window.end(), { X( 0 ), X( 1 ), X( 2 ) } );
  }
  return window;
}

/// the analytic tangent (material, fade and geometric parts) against central differences of the forces, with a
/// softening law with Coulomb friction and a window of 3 quadratic bar elements in a hexahedron, after a slide that
/// bends the bar and leaves the first channel point partially covered (0 < f < 1)
void testLargeSlipTangent3D()
{
  const Eigen::Vector3d a( 0.3, 0.2, 0.4 ), b( 1.6, 1.7, 1.5 );
  Bond                  large( "EBLS3D3H8W3",
              concat( hexa8, straightWindow3D( a, b, 7 ) ),
                               { 1.5, 0.0, 1.0, 1.0, 0.0, 4 },
              "MODELCODE2010BONDSLIP",
              mc2010Friction );
  const Eigen::Vector3d t    = ( b - a ).normalized();
  const double          L    = ( b - a ).norm() / 3.0; // bar element length = channel length
  Eigen::VectorXd       UOld = Eigen::VectorXd::Zero( large.nDof );
  for ( int i = 1; i <= 3; i++ ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( large.nDof );
    for ( int n = 0; n < 7; n++ )
      U.segment< 3 >( 24 + 3 * n ) = t * 0.1 * L * i;
    large.kernels( U, UOld, true );
    UOld = U;
  }
  Eigen::VectorXd U = UOld;
  // bent bar (curved quadratic elements), deformed host
  for ( int n = 0; n < 7; n++ )
    U.segment< 3 >( 24 + 3 * n ) += t * 0.05 * L + Eigen::Vector3d( 0.002, -0.001, 0.0005 ) * ( n * ( 6 - n ) );
  for ( int j = 0; j < 24; j++ )
    U( j ) += 0.001 * std::sin( 1.0 + j );

  const double covered = *large.element->getStateView( "covered fraction", 1 ).stateLocation;
  throwExceptionOnFailure( covered > 0 && covered < 1, "the second channel point must be partially covered" );

  const auto [P, K] = large.kernels( U, UOld, false );
  const auto numK   = large.numericalTangent( U, UOld, 1e-7 );
  throwExceptionOnFailure( checkIfEqual< double >( K, numK, 1e-5 * numK.cwiseAbs().maxCoeff() ),
                           Marmot::Utilities::MakeString()
                             << "large-slip tangent does not match central differences, max difference "
                             << ( K - numK ).cwiseAbs().maxCoeff() << " of " << numK.cwiseAbs().maxCoeff() );
}

/// a rigid rotation of the deformed configuration (host and bar) rotates the forces: P(Q x) = Q P(x), in 2D with
/// the MC2010 law with friction and in 3D with a linear law (LS-01)
void testLargeSlipObjectivity()
{
  {
    const std::vector< double > window = { 0.5, 1.5, 2.0, 1.6, 3.5, 1.5 };
    Bond                        large( "EBLS2D2Q8W2",
                concat( quad8, window ),
                                       { 1.0, 0.2, 0.9, 1.0, 1.0 },
                "MODELCODE2010BONDSLIP",
                mc2010Friction );
    const auto                  X = Eigen::Map< const Eigen::VectorXd >( large.coordinates.data(), large.nDof );
    Eigen::VectorXd             U = Eigen::VectorXd::Zero( large.nDof );
    for ( int n = 0; n < 3; n++ )
      U.segment< 2 >( 16 + 2 * n ) = Eigen::Vector2d( 0.3, 0.02 * n );
    U.head( 16 ) = Eigen::VectorXd::LinSpaced( 16, 0.0, 0.01 );
    for ( const double angle : { 0.3, 1.5708, 2.5 } ) {
      const Eigen::Rotation2Dd Q( angle );
      Eigen::VectorXd          URotated( large.nDof );
      for ( int n = 0; n < large.nDof / 2; n++ )
        URotated.segment< 2 >( 2 * n ) = Q * ( X.segment< 2 >( 2 * n ) + U.segment< 2 >( 2 * n ) ) -
                                         X.segment< 2 >( 2 * n );
      const Eigen::VectorXd P        = large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), false ).first;
      const Eigen::VectorXd PRotated = large.kernels( URotated, Eigen::VectorXd::Zero( large.nDof ), false ).first;
      Eigen::VectorXd       QP( large.nDof );
      for ( int n = 0; n < large.nDof / 2; n++ )
        QP.segment< 2 >( 2 * n ) = Q * P.segment< 2 >( 2 * n );
      throwExceptionOnFailure( P.norm() > 1 && checkIfEqual< double >( PRotated, QP, 1e-9 * P.norm() ),
                               Marmot::Utilities::MakeString()
                                 << "2D: P(Q x) != Q P(x) at " << angle
                                 << ", |difference| / |P| = " << ( PRotated - QP ).norm() / P.norm() );
    }
  }
  {
    const Eigen::Vector3d a( 0.3, 0.2, 0.4 ), b( 1.6, 1.7, 1.5 );
    Bond                  large( "EBLS3D2H8W3",
                concat( hexa8, straightWindow3D( a, b, 4 ) ),
                                 { 1.5, 1.0, 2.0, 1.0, 1.0 },
                "LINEARELASTICBONDSLIP",
                linearBond );
    const auto            X = Eigen::Map< const Eigen::VectorXd >( large.coordinates.data(), large.nDof );
    Eigen::VectorXd       U = Eigen::VectorXd::Zero( large.nDof );
    for ( int n = 0; n < 4; n++ )
      U.segment< 3 >( 24 + 3 * n ) = ( b - a ).normalized() * 0.2 + Eigen::Vector3d( 0.01, -0.02, 0.005 );
    const Eigen::Matrix3d Q = Eigen::AngleAxisd( 1.2, Eigen::Vector3d( 1, 2, -1 ).normalized() ).toRotationMatrix();
    Eigen::VectorXd       URotated( large.nDof );
    for ( int n = 0; n < large.nDof / 3; n++ )
      URotated.segment< 3 >( 3 * n ) = Q * ( X.segment< 3 >( 3 * n ) + U.segment< 3 >( 3 * n ) ) -
                                       X.segment< 3 >( 3 * n );
    const Eigen::VectorXd P        = large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), false ).first;
    const Eigen::VectorXd PRotated = large.kernels( URotated, Eigen::VectorXd::Zero( large.nDof ), false ).first;
    Eigen::VectorXd       QP( large.nDof );
    for ( int n = 0; n < large.nDof / 3; n++ )
      QP.segment< 3 >( 3 * n ) = Q * P.segment< 3 >( 3 * n );
    throwExceptionOnFailure( P.norm() > 1 && checkIfEqual< double >( PRotated, QP, 1e-9 * P.norm() ),
                             Marmot::Utilities::MakeString() << "3D: P(Q x) != Q P(x), |difference| / |P| = "
                                                             << ( PRotated - QP ).norm() / P.norm() );
  }
}

/// a host under uniform stretch along the bar and a rigid rotation, with the bar moving with it: no slip, no force
void testLargeSlipHostStretch()
{
  const std::vector< double > window = { 0.5, 1.5, 2.0, 1.5, 3.5, 1.5 };
  Bond                        large( "EBLS2D2Q8W2",
              concat( quad8, window ),
                                     { 1.0, 0.2, 0.9, 1.0, 1.0 },
              "LINEARELASTICBONDSLIP",
              linearBond );
  const auto                  X = Eigen::Map< const Eigen::VectorXd >( large.coordinates.data(), large.nDof );
  Eigen::Matrix2d             F;
  F << 1.2, 0.0, 0.0, 0.9; // stretch along the bar (x), contraction across
  F = Eigen::Rotation2Dd( 0.4 ).toRotationMatrix() * F;
  Eigen::VectorXd U( large.nDof );
  for ( int n = 0; n < large.nDof / 2; n++ )
    U.segment< 2 >( 2 * n ) = ( F - Eigen::Matrix2d::Identity() ) * X.segment< 2 >( 2 * n );
  const auto P = large.kernels( U, Eigen::VectorXd::Zero( large.nDof ), true ).first;
  // round-off of the closest point iteration only: a slip of 1e-9 times the displacements would give Kn 1e-9 |U|
  throwExceptionOnFailure( P.cwiseAbs().maxCoeff() < 1e-9 * 1000. * U.cwiseAbs().maxCoeff(),
                           Marmot::Utilities::MakeString() << "uniform stretch must not slip, |P| = " << P.norm() );
}

/// a bar pulled out until it leaves the first channel point (inactive, history frozen) and pushed back: the point is
/// active again and its law continues from the frozen history (MC2010: the largest slip)
void testLargeSlipPushOutPushBack()
{
  const std::vector< double > bar = { 0.5, 1.5, 3.5, 1.5 };
  Bond large( "EBLS2D2Q8W1", concat( quad8, bar ), { 1.0, 0.0, 1.0, 1.0, 1.0, 3 }, "MODELCODE2010BONDSLIP", mc2010 );
  auto maxSlip = [&]( int point ) { return *large.element->getStateView( "max slip", point ).stateLocation; };
  auto active  = [&]( int point ) { return *large.element->getStateView( "active", point ).stateLocation; };

  Eigen::VectorXd UOld      = Eigen::VectorXd::Zero( large.nDof );
  auto            moveBarTo = [&]( double u ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( large.nDof );
    U( 16 ) = U( 18 ) = u;
    large.kernels( U, UOld, true );
    UOld = U;
  };
  for ( int i = 1; i <= 12; i++ )
    moveBarTo( 0.05 * i ); // 0.6: the bar start passed the first channel point (x = 0.5)
  throwExceptionOnFailure( active( 0 ) == 0.0, "the first channel point must be inactive after the pull-out" );
  const double frozen = maxSlip( 0 );
  throwExceptionOnFailure( frozen > 0.3 && frozen < 0.6,
                           Marmot::Utilities::MakeString() << "frozen history " << frozen );
  moveBarTo( 0.7 );
  throwExceptionOnFailure( maxSlip( 0 ) == frozen, "the history of an inactive point must be frozen" );
  for ( int i = 1; i <= 12; i++ )
    moveBarTo( 0.7 - 0.05 * i ); // back to 0.1
  throwExceptionOnFailure( active( 0 ) == 1.0, "the first channel point must be active after the push-back" );
  throwExceptionOnFailure( maxSlip( 0 ) == frozen, "the law must continue from the frozen largest slip" );
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
          "EB2D2T3",
          "EB2D3T3",
          "EB2D2Q4",
          "EB2D2Q8",
          "EB2D3Q4",
          "EB2D3Q8",
          "EB3D2T4",
          "EB3D2H8",
          "EB3D2H20",
          "EB3D3T4",
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
    testOffsetCoordinates,
    testTriangleHost,
    testInvalidElementProperties,
    testTangentWithFriction,
    testLargeSlipEqualsSmallSlipForSmallDisplacements,
    testLargeSlipPullOutLosesBond,
    testLargeSlipWindowExceededThrows,
    testLargeSlipTangent3D,
    testLargeSlipUnequalBarElements,
    testLargeSlipObjectivity,
    testLargeSlipHostStretch,
    testLargeSlipPushOutPushBack,
  } );
  return 0;
}
