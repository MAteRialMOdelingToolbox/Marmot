#include "Marmot/BeamElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Discretizations::ElementProperties;

namespace {

  /// A single beam element with its own state, for tests
  struct Beam {
    std::unique_ptr< Marmot::Discretizations::Elements::Element > element;
    std::vector< double >                                         coordinates;
    std::vector< double >                                         properties;
    std::vector< double >                                         materialProperties;
    std::vector< double >                                         stateVars;
    int                                                           nDof;

    Beam( const std::string&    type,
          std::vector< double > coordinates_,
          std::vector< double > properties_,
          const std::string&    material,
          std::vector< double > materialProperties_ )
      : element( Marmot::Registry::ElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( std::move( properties_ ) ),
        materialProperties( std::move( materialProperties_ ) )
    {
      nDof = element->getNDofPerElement();
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty( ElementProperties( properties.data(), static_cast< int >( properties.size() ) ) );
      element->initializeYourself();
      element->assignProperty(
        Marmot::Discretizations::MaterialSection( material,
                                                  materialProperties.data(),
                                                  static_cast< int >( materialProperties.size() ) ) );
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->setInitialConditions( Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization,
                                     nullptr );
    }

    /// residual and tangent for the total dofs U from the committed state; commits if requested
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

    /// central difference tangent from the committed state
    Eigen::MatrixXd numericalTangent( const Eigen::VectorXd& U, const Eigen::VectorXd& UOld, double h )
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

    /// the linear elastic solution of the first node clamped, the second node (the other end) loaded by f, plus the
    /// nodal loads fAll (all nodes, e.g., of a body force); returns the dofs of the second node
    Eigen::VectorXd cantilever( const Eigen::VectorXd& f, const Eigen::VectorXd& fAll = Eigen::VectorXd() )
    {
      const int       n = nDof / element->getNNodes(); // dofs per node
      const auto      K = kernels( Eigen::VectorXd::Zero( nDof ), Eigen::VectorXd::Zero( nDof ), false ).second;
      Eigen::VectorXd F = fAll.size() ? Eigen::VectorXd( fAll.tail( nDof - n ) ) : Eigen::VectorXd::Zero( nDof - n );
      F.head( n ) += f;
      const Eigen::VectorXd u = K.bottomRightCorner( nDof - n, nDof - n ).ldlt().solve( F );
      return u.head( n );
    }

    Eigen::VectorXd state( const std::string& name, int qp = 0 )
    {
      const auto view = element->getStateView( name, qp );
      return Eigen::Map< const Eigen::VectorXd >( view.stateLocation, view.stateSize );
    }
  };

  const double                E = 210000., nu = 0.3, G = E / ( 2 * ( 1 + nu ) );
  const std::vector< double > linearElastic = { E, nu, 7.85e-9 };

  void check( double value, double expected, double relTol, const std::string& what )
  {
    throwExceptionOnFailure( std::abs( value - expected ) <= relTol * std::abs( expected ),
                             Marmot::Utilities::MakeString() << what << ": " << value << " != " << expected );
  }

  /// drive the beam from U = 0 to U in nSteps equal, accepted increments, as an analysis would: the material sees
  /// moderate trial stresses (the return mapping of a nearly perfectly plastic material is then well conditioned on
  /// every platform)
  void loadInSteps( Beam& beam, const Eigen::VectorXd& U, int nSteps )
  {
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( U.size() );
    for ( int i = 1; i <= nSteps; i++ ) {
      const Eigen::VectorXd Ui = U * ( double( i ) / nSteps );
      beam.kernels( Ui, UOld, true );
      UOld = Ui;
    }
  }

  void checkTangent( Beam&                  beam,
                     const Eigen::VectorXd& U,
                     const Eigen::VectorXd& UOld,
                     double                 h,
                     double                 relTol,
                     const std::string&     what )
  {
    const auto [P, K]  = beam.kernels( U, UOld, false );
    const auto   numK  = beam.numericalTangent( U, UOld, h );
    const double scale = numK.cwiseAbs().maxCoeff();
    throwExceptionOnFailure( checkIfEqual< double >( K, numK, relTol * scale ),
                             Marmot::Utilities::MakeString()
                               << what << ": tangent does not match the numerical one\nK=\n"
                               << K << "\nnumK=\n"
                               << numK );
  }

  /// section points (y, z, A) of a rectangle y in [y0, y1], z in [z0, z1]: ny x nz cells with 2 x 2 Gauss points each
  /// (exact area, first and second moments)
  std::vector< std::array< double, 3 > > rectanglePoints( double y0, double y1, double z0, double z1, int ny, int nz )
  {
    const double                           g = 1 / std::sqrt( 3. );
    std::vector< std::array< double, 3 > > points;
    const double                           dy = ( y1 - y0 ) / ny, dz = ( z1 - z0 ) / nz;
    for ( int i = 0; i < ny; i++ )
      for ( int j = 0; j < nz; j++ )
        for ( double a : { -g, g } )
          for ( double b : { -g, g } )
            points.push_back( { y0 + ( i + 0.5 + 0.5 * a ) * dy, z0 + ( j + 0.5 + 0.5 * b ) * dz, 0.25 * dy * dz } );
    return points;
  }

  /// section points of a circle (diameter d): rings with 2 radial Gauss points (weight r) x uniform sectors (exact A,
  /// first and second moments)
  std::vector< std::array< double, 3 > > circlePoints( double d, int nRings, int nSectors )
  {
    const double                           g = 1 / std::sqrt( 3. ), R = d / 2;
    std::vector< std::array< double, 3 > > points;
    for ( int r = 0; r < nRings; r++ ) {
      const double r0 = R * r / nRings, r1 = R * ( r + 1 ) / nRings;
      for ( double a : { -g, g } ) {
        const double rho = 0.5 * ( r0 + r1 ) + 0.5 * a * ( r1 - r0 );
        for ( int k = 0; k < nSectors; k++ ) {
          const double phi = ( k + 0.5 ) * 2 * std::numbers::pi / nSectors;
          points.push_back( { rho * std::cos( phi ),
                              rho * std::sin( phi ),
                              rho * 0.5 * ( r1 - r0 ) * 2 * std::numbers::pi / nSectors } );
        }
      }
    }
    return points;
  }

  /// the element properties of a B23 for section points (z dropped)
  std::vector< double > properties2D( const std::vector< std::array< double, 3 > >& points )
  {
    std::vector< double > p = { double( points.size() ) };
    for ( const auto& q : points )
      p.insert( p.end(), { q[0], q[2] } );
    return p;
  }

  /// the element properties of a B33 for section points, orientation and torsion constant
  std::vector< double > properties3D( const std::vector< std::array< double, 3 > >& points,
                                      const Eigen::Vector3d&                        v,
                                      double                                        J )
  {
    std::vector< double > p = { v( 0 ), v( 1 ), v( 2 ), J, double( points.size() ) };
    for ( const auto& q : points )
      p.insert( p.end(), { q[0], q[1], q[2] } );
    return p;
  }

  /// nodal coordinates of a straight beam from X1 to X2 with 2 or 3 nodes (end, end, mid)
  std::vector< double > beamNodes( const std::vector< double >& X1, const std::vector< double >& X2, int nNodes )
  {
    std::vector< double > X = X1;
    X.insert( X.end(), X2.begin(), X2.end() );
    if ( nNodes == 3 )
      for ( size_t i = 0; i < X1.size(); i++ )
        X.push_back( 0.5 * ( X1[i] + X2[i] ) );
    return X;
  }

  std::vector< std::array< double, 3 > > shifted( std::vector< std::array< double, 3 > > points, double dy, double dz )
  {
    for ( auto& q : points )
      q[0] += dy, q[1] += dz;
    return points;
  }

  std::vector< std::array< double, 3 > > concat( std::vector< std::array< double, 3 > >        a,
                                                 const std::vector< std::array< double, 3 > >& b )
  {
    a.insert( a.end(), b.begin(), b.end() );
    return a;
  }

  /// analytical integrals of a union of rectangles {y0, y1, z0, z1}: A, int y, int z, int z^2, int y^2, int y z
  std::array< double, 6 > rectanglesIntegrals( const std::vector< std::array< double, 4 > >& rectangles )
  {
    std::array< double, 6 > c{};
    for ( const auto& [y0, y1, z0, z1] : rectangles ) {
      const double Y1 = y1 - y0, Y2 = ( y1 * y1 - y0 * y0 ) / 2, Y3 = ( y1 * y1 * y1 - y0 * y0 * y0 ) / 3;
      const double Z1 = z1 - z0, Z2 = ( z1 * z1 - z0 * z0 ) / 2, Z3 = ( z1 * z1 * z1 - z0 * z0 * z0 ) / 3;
      c[0] += Y1 * Z1, c[1] += Y2 * Z1, c[2] += Y1 * Z2, c[3] += Y1 * Z3, c[4] += Y3 * Z1, c[5] += Y2 * Z2;
    }
    return c;
  }

} // namespace

/// the beams are registered, carry displacement and rotation, and report a bar2 shape
void testFactoryNames()
{
  for ( const auto& [name, nDim, nDofPerNode] :
        std::vector< std::tuple< std::string, int, int > >{ { "B2D2", 2, 3 }, { "B3D2", 3, 6 } } ) {
    std::unique_ptr< Marmot::Discretizations::Elements::Element > el(
      Marmot::Registry::ElementFactory::createElement( name, 1 ) );
    throwExceptionOnFailure( el->getNSpatialDimensions() == nDim && el->getNDofPerElement() == 2 * nDofPerNode &&
                               el->getElementShape() == "bar2" && el->getNNodes() == 2,
                             name + " element description" );
    const auto fields = el->getNodeFields();
    throwExceptionOnFailure( fields.size() == 2 &&
                               fields[0] == std::vector< std::string >{ "displacement", "rotation" },
                             name + " node fields" );
  }
}

/// 2D cantilever (inclined, a single element) under a tip force, a tip moment and an axial force: exact for a
/// rectangle and a circle (exactly integrated section points); the moment at the quadrature points is P (L - x)
void testCantilever2D()
{
  const double          b = 10, h = 20, d = 20, L = 500;
  const Eigen::Vector2d t( 0.6, 0.8 ), n( -0.8, 0.6 );
  const std::vector< std::tuple< std::string, std::vector< std::array< double, 3 > >, double, double > > sections =
    { { "rectangle", rectanglePoints( -h / 2, h / 2, -b / 2, b / 2, 3, 1 ), b * h, b * h * h * h / 12 },
      { "circle", circlePoints( d, 1, 16 ), std::numbers::pi * d * d / 4, std::numbers::pi * d * d * d * d / 64 } };
  for ( int nNodes : { 2, 3 } )
    for ( const auto& [name, points, A, I] : sections ) {
      const std::string type = "B2D" + std::to_string( nNodes );
      const auto        X    = beamNodes( { 100., 50. }, { 400., 450. }, nNodes );
      Beam              beam( type, X, properties2D( points ), "LINEARELASTIC", linearElastic );
      const std::string what = type + " " + name + ": ";
      const int         nDof = 3 * nNodes;

      const double    P = 7.5, M = 1200., N = 3000.;
      Eigen::VectorXd f( 3 );
      f.head< 2 >()     = P * n;
      f( 2 )            = 0;
      Eigen::VectorXd u = beam.cantilever( f );
      check( u.head( 2 ).dot( n ), P * L * L * L / ( 3 * E * I ), 1e-10, what + "tip deflection (force)" );
      check( u( 2 ), P * L * L / ( 2 * E * I ), 1e-10, what + "tip rotation (force)" );
      throwExceptionOnFailure( std::abs( u.head( 2 ).dot( t ) ) < 1e-12, what + "no axial displacement (force)" );

      // the full solution (the mid node of the 3-node beam from the exact cubic deflection)
      Eigen::VectorXd U = Eigen::VectorXd::Zero( nDof );
      U.segment( 3, 3 ) = u;
      if ( nNodes == 3 ) {
        const double vMid   = P * ( L / 2 ) * ( L / 2 ) * ( 3 * L - L / 2 ) / ( 6 * E * I );
        const double tMid   = P * ( L / 2 ) * ( 2 * L - L / 2 ) / ( 2 * E * I );
        U.segment< 2 >( 6 ) = vMid * n;
        U( 8 )              = tMid;
      }
      beam.kernels( U, Eigen::VectorXd::Zero( nDof ), true );
      const auto xs = beam.element->getCoordinatesAtQuadraturePoints();
      for ( int qp = 0; qp < beam.element->getNumberOfQuadraturePoints(); qp++ ) {
        const double x = ( Eigen::Vector2d( xs[qp][0], xs[qp][1] ) - Eigen::Vector2d( 100., 50. ) ).norm();
        throwExceptionOnFailure( std::abs( beam.state( "bending moment", qp )( 0 ) - P * ( L - x ) ) < 1e-9 * P * L,
                                 what + "bending moment" );
        throwExceptionOnFailure( std::abs( beam.state( "curvature", qp )( 0 ) - P * ( L - x ) / ( E * I ) ) <
                                   1e-9 * P * L / ( E * I ),
                                 what + "curvature" );
      }
      {
        // the linear elastic material does not track its energy; an elastic von Mises material does
        Beam vm( type, X, properties2D( points ), "VONMISES", { E, nu, 1e12, 0., 0., 1., 7.85e-9 } );
        vm.kernels( U, Eigen::VectorXd::Zero( nDof ), true );
        double energy = 0;
        vm.element->computeInternalEnergy( energy );
        check( energy, 0.5 * P * u.head( 2 ).dot( n ), 1e-10, what + "elastic energy" );
      }

      f << 0, 0, M;
      u = beam.cantilever( f );
      check( u.head( 2 ).dot( n ), M * L * L / ( 2 * E * I ), 1e-10, what + "tip deflection (moment)" );
      check( u( 2 ), M * L / ( E * I ), 1e-10, what + "tip rotation (moment)" );

      f.head< 2 >() = N * t;
      f( 2 )        = 0;
      u             = beam.cantilever( f );
      check( u.head( 2 ).dot( t ), N * L / ( E * A ), 1e-10, what + "axial elongation" );
      throwExceptionOnFailure( std::abs( u( 2 ) ) < 1e-14, what + "no rotation (axial)" );

      // a uniform line load q (body force b A): the 3-node beam (quintic) is exact, q L^4 / (8 E I)
      if ( nNodes == 3 ) {
        const double          q  = 0.01;
        Eigen::VectorXd       fb = Eigen::VectorXd::Zero( nDof );
        const Eigen::Vector2d bf = q / A * n;
        beam.element->computeBodyForce( fb.data(), nullptr, bf.data(), nullptr, 0, 1 );
        u = beam.cantilever( Eigen::VectorXd::Zero( 3 ), fb );
        check( u.head( 2 ).dot( n ), q * L * L * L * L / ( 8 * E * I ), 1e-10, what + "tip deflection (line load)" );
        check( u( 2 ), q * L * L * L / ( 6 * E * I ), 1e-10, what + "tip rotation (line load)" );
      }
    }
}

/// 3D cantilever in a general orientation: bending about both local axes and torsion (rectangle and circle)
void testCantilever3D()
{
  const double          L = 300, b = 10, h = 20, d = 20;
  const Eigen::Vector3d X1( 10., -20., 5. ), t = Eigen::Vector3d( 1., 2., 2. ) / 3.;
  const Eigen::Vector3d X2 = X1 + L * t;
  const Eigen::Vector3d v( 0., 0., 1. );
  const Eigen::Vector3d e2 = ( v - v.dot( t ) * t ).normalized(), e3 = t.cross( e2 );

  // name, points, Iy, Iz, J
  const std::vector< std::tuple< std::string, std::vector< std::array< double, 3 > >, double, double, double > >
    sections = { { "rectangle",
                   rectanglePoints( -h / 2, h / 2, -b / 2, b / 2, 2, 2 ),
                   h * b * b * b / 12,
                   b * h * h * h / 12,
                   2500. },
                 { "circle",
                   circlePoints( d, 1, 16 ),
                   std::numbers::pi * d * d * d * d / 64,
                   std::numbers::pi * d * d * d * d / 64,
                   std::numbers::pi * d * d * d * d / 32 } };
  for ( int nNodes : { 2, 3 } )
    for ( const auto& [name, points, Iy, Iz, J] : sections ) {
      const std::string type = "B3D" + std::to_string( nNodes );
      Beam              beam( type,
                 beamNodes( { X1( 0 ), X1( 1 ), X1( 2 ) }, { X2( 0 ), X2( 1 ), X2( 2 ) }, nNodes ),
                 properties3D( points, v, J ),
                 "LINEARELASTIC",
                 linearElastic );
      const std::string what = type + " " + name + ": ";
      const double      P = 5., Tq = 2000.;

      Eigen::VectorXd f( 6 ), u;
      f.head< 3 >() = P * e2;
      f.tail< 3 >().setZero();
      u = beam.cantilever( f );
      check( u.head( 3 ).dot( e2 ), P * L * L * L / ( 3 * E * Iz ), 1e-10, what + "deflection along e2" );
      check( u.tail( 3 ).dot( e3 ), P * L * L / ( 2 * E * Iz ), 1e-10, what + "rotation about e3" );
      throwExceptionOnFailure( std::abs( u.head( 3 ).dot( e3 ) ) + std::abs( u.tail( 3 ).dot( e2 ) ) < 1e-12,
                               what + "no coupling e2 -> e3" );

      f.head< 3 >() = P * e3;
      f.tail< 3 >().setZero();
      u = beam.cantilever( f );
      check( u.head( 3 ).dot( e3 ), P * L * L * L / ( 3 * E * Iy ), 1e-10, what + "deflection along e3" );
      check( u.tail( 3 ).dot( e2 ), -P * L * L / ( 2 * E * Iy ), 1e-10, what + "rotation about e2" );

      f.head< 3 >().setZero();
      f.tail< 3 >() = Tq * t;
      u             = beam.cantilever( f );
      check( u.tail( 3 ).dot( t ), Tq * L / ( G * J ), 1e-10, what + "twist" );
      throwExceptionOnFailure( u.head( 3 ).norm() < 1e-12, what + "no displacement under torsion" );

      // torsion: the twist is linear along the beam
      Eigen::VectorXd U = Eigen::VectorXd::Zero( 6 * nNodes );
      U.segment( 6, 6 ) = u;
      if ( nNodes == 3 )
        U.segment( 12, 6 ) = 0.5 * u;
      beam.kernels( U, Eigen::VectorXd::Zero( 6 * nNodes ), true );
      check( beam.state( "torque", 1 )( 0 ), Tq, 1e-10, what + "torque" );
      check( beam.state( "twist rate", 1 )( 0 ), Tq / ( G * J ), 1e-10, what + "twist rate" );
    }
}

/// a non-symmetric L-profile (centroid on the axis): a tip force along the local y axis, which is not a principal
/// axis, also deflects the beam along z; exact with the analytical section constants Iy, Iz, Iyz
void testUnsymmetricLProfile3D()
{
  const double a = 60, bLeg = 40, th = 6, L = 400;
  // legs: y in [0, a] x z in [0, th], and y in [0, th] x z in [th, bLeg]
  const std::vector< std::array< double, 4 > > rectangles = { { 0, a, 0, th }, { 0, th, th, bLeg } };
  const auto                                   c          = rectanglesIntegrals( rectangles );
  const double                                 A = c[0], yc = c[1] / A, zc = c[2] / A;
  const double Iy = c[3] - A * zc * zc, Iz = c[4] - A * yc * yc, Iyz = c[5] - A * yc * zc;
  throwExceptionOnFailure( std::abs( Iyz ) > 0.1 * std::sqrt( Iy * Iz ), "the L-profile axes must not be principal" );

  const auto points = shifted( concat( rectanglePoints( 0, a, 0, th, 6, 1 ), rectanglePoints( 0, th, th, bLeg, 1, 4 ) ),
                               -yc,
                               -zc );
  Beam       beam( "B3D2",
                   { 0., 0., 0., L, 0., 0. },
             properties3D( points, { 0., 1., 0. }, A * th * th / 3 ),
             "LINEARELASTIC",
             linearElastic );

  // the integrated section constants of the points
  const auto* section = &static_cast< Marmot::Discretizations::Elements::BeamElement< 3, 2 >* >( beam.element.get() )
                           ->qps[0]
                           .section;
  const auto ci = section->integrals();
  check( ci[0], A, 1e-12, "L-profile: A of the points" );
  check( ci[3], Iy, 1e-12, "L-profile: Iy of the points" );
  check( ci[4], Iz, 1e-12, "L-profile: Iz of the points" );
  check( ci[5], Iyz, 1e-12, "L-profile: Iyz of the points" );

  const double    Fy = 3., Fz = -1.;
  Eigen::VectorXd f( 6 );
  f << 0, Fy, Fz, 0, 0, 0;
  const Eigen::VectorXd u = beam.cantilever( f );
  // [My, Mz] = E [[Iy, -Iyz], [-Iyz, Iz]] [kappa_y, kappa_z], My = -Fz (L - x), Mz = Fy (L - x)
  Eigen::Matrix2d Db;
  Db << Iy, -Iyz, -Iyz, Iz;
  const Eigen::Vector2d kappa = ( E * Db ).inverse() * Eigen::Vector2d( -Fz, Fy ); // per (L - x)
  check( u( 1 ), kappa( 1 ) * L * L * L / 3, 1e-9, "L-profile: deflection along y" );
  check( u( 2 ), -kappa( 0 ) * L * L * L / 3, 1e-9, "L-profile: deflection along z (coupled)" );
  throwExceptionOnFailure( std::abs( u( 0 ) ) < 1e-12, "L-profile: no axial displacement (centroid on the axis)" );
}

/// (infinitesimal) rigid body motions are stress free
void testRigidBodyMotions()
{
  const auto points = rectanglePoints( -5, 5, -2, 2, 2, 2 );
  {
    Beam            beam( "B2D2", { 1., 2., 41., 32. }, properties2D( points ), "LINEARELASTIC", linearElastic );
    const double    w = 1e-3;
    Eigen::VectorXd U( 6 );
    U << 0.3 - w * 2., 0.1 + w * 1., w, 0.3 - w * 32., 0.1 + w * 41., w;
    const auto P = beam.kernels( U, Eigen::VectorXd::Zero( 6 ), false ).first;
    throwExceptionOnFailure( P.norm() < 1e-10,
                             Marmot::Utilities::MakeString() << "B23 rigid body motion: |P| = " << P.norm() );
  }
  {
    Beam                  beam( "B3D2",
                                { 1., 2., 3., 41., 32., -7. },
               properties3D( points, { 0., 0., 1. }, 80. ),
               "LINEARELASTIC",
               linearElastic );
    const Eigen::Vector3d c( 0.3, -0.2, 0.1 ), w( 1e-3, -2e-3, 0.5e-3 );
    Eigen::VectorXd       U( 12 );
    U.segment< 3 >( 0 ) = c + w.cross( Eigen::Vector3d( 1., 2., 3. ) );
    U.segment< 3 >( 3 ) = w;
    U.segment< 3 >( 6 ) = c + w.cross( Eigen::Vector3d( 41., 32., -7. ) );
    U.segment< 3 >( 9 ) = w;
    const auto P        = beam.kernels( U, Eigen::VectorXd::Zero( 12 ), false ).first;
    throwExceptionOnFailure( P.norm() < 1e-10,
                             Marmot::Utilities::MakeString() << "B33 rigid body motion: |P| = " << P.norm() );
  }
  {
    Beam                  beam( "B3D3",
               beamNodes( { 1., 2., 3. }, { 41., 32., -7. }, 3 ),
               properties3D( points, { 0., 0., 1. }, 80. ),
               "LINEARELASTIC",
               linearElastic );
    const Eigen::Vector3d c( 0.3, -0.2, 0.1 ), w( 1e-3, -2e-3, 0.5e-3 );
    Eigen::VectorXd       U( 18 );
    U.segment< 3 >( 0 )  = c + w.cross( Eigen::Vector3d( 1., 2., 3. ) );
    U.segment< 3 >( 3 )  = w;
    U.segment< 3 >( 6 )  = c + w.cross( Eigen::Vector3d( 41., 32., -7. ) );
    U.segment< 3 >( 9 )  = w;
    U.segment< 3 >( 12 ) = c + w.cross( Eigen::Vector3d( 21., 17., -2. ) );
    U.segment< 3 >( 15 ) = w;
    const auto P         = beam.kernels( U, Eigen::VectorXd::Zero( 18 ), false ).first;
    throwExceptionOnFailure( P.norm() < 1e-8,
                             Marmot::Utilities::MakeString() << "B3D3 rigid body motion: |P| = " << P.norm() );
  }
}

/// consistent nodal loads of a uniform line load: total force b A L and fixed-end moments q L^2 / 12
void testBodyForce()
{
  const double    A = 20, L = 50;
  Beam            beam( "B2D2",
                        { 0., 0., 30., 40. },
             properties2D( rectanglePoints( -2, 2, -2.5, 2.5, 1, 1 ) ),
             "LINEARELASTIC",
             linearElastic );
  Eigen::VectorXd P    = Eigen::VectorXd::Zero( 6 );
  const double    b[2] = { 0., -2. };
  beam.element->computeBodyForce( P.data(), nullptr, b, nullptr, 0, 1 );
  check( P( 1 ) + P( 4 ), -2. * A * L, 1e-12, "body force: total" );
  const double q = -2. * A * 0.6; // transverse component (normal (-0.8, 0.6))
  check( P( 2 ), q * L * L / 12, 1e-12, "body force: moment at node 1" );
  check( P( 5 ), -q * L * L / 12, 1e-12, "body force: moment at node 2" );
}

/// fiber section with von Mises plasticity: consistent tangents in plastic states (bending, axial force, torsion),
/// also for the unsymmetric L-profile
void testPlasticTangents()
{
  //                                     E        nu   fy    H     dfy  delta rho
  const std::vector< double > vm = { 210000., 0.3, 500., 2000., 0., 1., 7.85e-9 };
  {
    Beam beam( "B2D2", { 0., 0., 60., 80. }, properties2D( rectanglePoints( -10, 10, -5, 5, 3, 1 ) ), "VONMISES", vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 6 ), U( 6 );
    U << 0., 0., 0.02, 0.15, 0.05, 0.06;
    loadInSteps( beam, U, 20 );
    UOld = U;
    U << 0., 0., 0.025, 0.17, 0.06, 0.07;
    checkTangent( beam, U, UOld, 1e-8, 1e-5, "B23 von Mises" );
  }
  const auto rectangle = rectanglePoints( -10, 10, -5, 5, 2, 2 );
  const auto lProfile  = concat( rectanglePoints( -20, 20, -3, 3, 4, 1 ), rectanglePoints( -20, -14, 3, 20, 1, 3 ) );
  for ( const auto& points : { rectangle, lProfile } ) {
    Beam beam( "B3D2", { 0., 0., 0., 60., 80., 0. }, properties3D( points, { 0., 0., 1. }, 4500. ), "VONMISES", vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 12 ), U( 12 );
    U << 0., 0., 0., 0., 0., 0., 0.15, 0.05, 0.1, 0.05, -0.02, 0.06;
    loadInSteps( beam, U, 20 );
    UOld = U;
    U << 0., 0., 0., 0., 0., 0., 0.17, 0.06, 0.12, 0.06, -0.025, 0.07;
    checkTangent( beam, U, UOld, 1e-8, 1e-5, "B33 von Mises (bending, axial, torsion)" );
    throwExceptionOnFailure( beam.state( "dissipation", 0 )( 0 ) > 0, "B33 von Mises: plastic dissipation" );
  }
  {
    Beam            beam( "B3D3",
               beamNodes( { 0., 0., 0. }, { 60., 80., 0. }, 3 ),
               properties3D( lProfile, { 0., 0., 1. }, 4500. ),
               "VONMISES",
               vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 18 ), U( 18 );
    U << 0., 0., 0., 0., 0., 0., 0.15, 0.05, 0.1, 0.005, -0.002, 0.006, 0.07, 0.04, 0.06, 0.003, -0.001, 0.003;
    loadInSteps( beam, U, 20 );
    UOld = U;
    U( 7 ) += 0.02, U( 11 ) += 0.001, U( 15 ) += 0.0005;
    checkTangent( beam, U, UOld, 1e-8, 1e-5, "B3D3 von Mises" );
  }
}

/// fully plastic section points (perfect von Mises plasticity): the plastic moment f_y sum |y_f| A_f of the points
/// under pure bending, and the axial yield force f_y A
void testPlasticLimits()
{
  const double L = 100, fy = 500;
  const auto   points = rectanglePoints( -10, 10, -5, 5, 10, 1 );
  double       Mp = 0, A = 0;
  for ( const auto& q : points )
    Mp += fy * std::abs( q[0] ) * q[2], A += q[2];
  const std::vector< double > vm = { E, nu, fy, 0., 0., 1., 7.85e-9 };

  Beam            beam( "B2D2", { 0., 0., L, 0. }, properties2D( points ), "VONMISES", vm );
  const double    kappaMax = 200 * 2 * fy / ( E * 20 ); // all points yielded
  Eigen::VectorXd UOld     = Eigen::VectorXd::Zero( 6 );
  // small steps (outer fibre strain increment 0.5 eps_y): the return mapping of a perfectly plastic material stays
  // well conditioned on every platform
  for ( int i = 1; i <= 100; i++ ) {
    const double    Theta = 0.5 * i / 100. * kappaMax * L;
    Eigen::VectorXd U( 6 );
    U << 0, 0, -Theta, 0, 0, Theta;
    beam.kernels( U, UOld, true );
    UOld = U;
  }
  for ( int qp = 0; qp < 3; qp++ ) {
    check( beam.state( "bending moment", qp )( 0 ), Mp, 1e-6, "plastic moment" );
    throwExceptionOnFailure( std::abs( beam.state( "normal force", qp )( 0 ) ) < 1e-6 * fy * A, "pure bending" );
  }

  Beam bar( "B2D2", { 0., 0., L, 0. }, properties2D( points ), "VONMISES", vm );
  UOld = Eigen::VectorXd::Zero( 6 );
  for ( int i = 1; i <= 10; i++ ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( 6 );
    U( 3 )            = i * 0.05;
    bar.kernels( U, UOld, true );
    UOld = U;
  }
  check( bar.state( "normal force", 1 )( 0 ), fy * A, 1e-8, "axial yield force" );
}

/// limit load of a 2D cantilever (tip load) of an elastic-perfectly plastic (von Mises) beam, a chain of B23 elements
/// under tip displacement control (axial displacement free), by Newton iterations on the assembled system
double cantileverLimitLoad( int                                           nElements,
                            const std::vector< std::array< double, 3 > >& points,
                            double                                        L,
                            double                                        fy,
                            double                                        vMax )
{
  const std::vector< double > vm = { E, nu, fy, 0., 0., 1., 7.85e-9 };
  std::vector< Beam >         beams;
  for ( int e = 0; e < nElements; e++ )
    beams.emplace_back( "B2D2",
                        std::vector< double >{ e * L / nElements, 0., ( e + 1 ) * L / nElements, 0. },
                        properties2D( points ),
                        "VONMISES",
                        vm );
  const int       nDof = 3 * ( nElements + 1 );
  Eigen::VectorXd U = Eigen::VectorXd::Zero( nDof ), UOld = U;
  const int       tip = 3 * nElements + 1; // tip deflection dof
  // free dofs: all but the clamped root (0, 1, 2) and the controlled tip deflection
  std::vector< int > free;
  for ( int i = 3; i < nDof; i++ )
    if ( i != tip )
      free.push_back( i );

  // adaptive load stepping: an increment whose Newton iteration fails (diverges, or a material requests a cutback)
  // is repeated with half the size
  double force = 0, convergedForce = 0, v = 0, dv = vMax / 20;
  while ( v < vMax * ( 1 - 1e-12 ) ) {
    const double vNew      = std::min( vMax, v + dv );
    bool         converged = false;
    try {
      // predictor: the free dofs follow the tip with the tangent of the last converged state
      {
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
        for ( int e = 0; e < nElements; e++ )
          K.block( 3 * e,
                   3 * e,
                   6,
                   6 ) += beams[e].kernels( UOld.segment( 3 * e, 6 ), UOld.segment( 3 * e, 6 ), false ).second;
        Eigen::MatrixXd Kff( free.size(), free.size() );
        Eigen::VectorXd Kft( free.size() );
        for ( size_t i = 0; i < free.size(); i++ ) {
          Kft( i ) = K( free[i], tip );
          for ( size_t j = 0; j < free.size(); j++ )
            Kff( i, j ) = K( free[i], free[j] );
        }
        const Eigen::VectorXd du = Kff.lu().solve( -Kft * ( vNew - v ) );
        for ( size_t i = 0; i < free.size(); i++ )
          U( free[i] ) = UOld( free[i] ) + du( i );
        U( tip ) = vNew;
      }
      for ( int it = 0; it < 25; it++ ) {
        Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
        for ( int e = 0; e < nElements; e++ ) {
          const auto [Pe, Ke] = beams[e].kernels( U.segment( 3 * e, 6 ), UOld.segment( 3 * e, 6 ), false );
          P.segment( 3 * e, 6 ) += Pe;
          K.block( 3 * e, 3 * e, 6, 6 ) += Ke;
        }
        Eigen::VectorXd r( free.size() );
        Eigen::MatrixXd Kff( free.size(), free.size() );
        for ( size_t i = 0; i < free.size(); i++ ) {
          r( i ) = P( free[i] );
          for ( size_t j = 0; j < free.size(); j++ )
            Kff( i, j ) = K( free[i], free[j] );
        }
        force = P( tip );
        if ( r.norm() < 1e-9 * std::max( 1.0, std::abs( force ) ) ) {
          converged = true;
          break;
        }
        const Eigen::VectorXd du = Kff.lu().solve( -r );
        for ( size_t i = 0; i < free.size(); i++ )
          U( free[i] ) += du( i );
      }
    }
    catch ( const Marmot::Utilities::StressUpdateFailed& ) {
    }
    if ( !converged ) {
      U  = UOld;
      dv = 0.5 * dv;
      if ( dv < vMax / 320 )
        break; // the plateau: the hinge is a mechanism, the tangent singular
      continue;
    }
    for ( int e = 0; e < nElements; e++ )
      beams[e].kernels( U.segment( 3 * e, 6 ), UOld.segment( 3 * e, 6 ), true );
    UOld           = U;
    v              = vNew;
    convergedForce = force;
  }
  return convergedForce;
}

/// plastic hinge of a von Mises cantilever: the limit load converges to M_p / L, M_p = f_y Z, with the number of
/// elements and section points; rectangle (Z = b h^2 / 4), circle (Z = d^3 / 6) and an unsymmetric T-section (the
/// plastic neutral axis is the equal-area axis, not the centroidal one: the free axial dof lets it shift)
void testPlasticHingeCantilever()
{
  const double L = 1000, fy = 235;

  // the limit load with nCoarse and 2 nCoarse elements; it converges linearly with the element size (the hinge forms at
  // the root quadrature point), so the extrapolation 2 P_fine - P_coarse must match M_p / L
  const int nCoarse = 4;
  auto      check =
    [&]( const std::string& name, const std::vector< std::array< double, 3 > >& points, double PLimit, double height ) {
      const double vMax = 10 * 2 * fy / ( E * height ) * L * L / 3;
      const double rC   = cantileverLimitLoad( nCoarse, points, L, fy, vMax ) / PLimit;
      const double rF   = cantileverLimitLoad( 2 * nCoarse, points, L, fy, vMax ) / PLimit;
      std::cout << name << ": P / (f_y Z / L) = " << rC << " (" << nCoarse << " elements), " << rF << " ("
                << 2 * nCoarse << " elements), extrapolated " << 2 * rF - rC << "\n";
      throwExceptionOnFailure( std::abs( rF - 1 ) < 0.6 * std::abs( rC - 1 ),
                               Marmot::Utilities::MakeString() << name << " plastic hinge: no convergence" );
      throwExceptionOnFailure( std::abs( 2 * rF - rC - 1 ) < 1e-2,
                               Marmot::Utilities::MakeString()
                                 << name << " plastic hinge: extrapolated P/P_lim = " << 2 * rF - rC );
    };

  // rectangle: the 2 x 2 Gauss points per cell give Z = b h^2 / 4 exactly
  const double b = 20, h = 40;
  check( "rectangle", rectanglePoints( -h / 2, h / 2, -b / 2, b / 2, 2, 1 ), fy * b * h * h / 4 / L, h );

  const double d = 40;
  check( "circle", circlePoints( d, 1, 16 ), fy * d * d * d / 6 / L, d );

  // T-section: flange y in [hw, hw + tf] x width bf, web y in [0, hw] x width tw; centroid on the axis. The plastic
  // neutral axis is the equal-area axis, not the centroidal one: the free axial dof lets it shift
  const double bf = 60, tf = 8, hw = 50, tw = 6;
  const auto   c = rectanglesIntegrals( { { hw, hw + tf, -bf / 2, bf / 2 }, { 0, hw, -tw / 2, tw / 2 } } );
  const double A = c[0], yc = c[1] / A;
  const double yp    = bf * tf >= A / 2 ? hw + tf - A / 2 / bf : ( A / 2 ) / tw;
  auto         Zpart = []( double y0, double y1, double w, double y ) {
    auto F = []( double s ) { return s * std::abs( s ) / 2; };
    return w * ( F( y1 - y ) - F( y0 - y ) );
  };
  const double Z = Zpart( hw, hw + tf, bf, yp ) + Zpart( 0, hw, tw, yp );
  check( "T-section",
         shifted( concat( rectanglePoints( hw, hw + tf, -bf / 2, bf / 2, 4, 1 ),
                          rectanglePoints( 0, hw, -tw / 2, tw / 2, 6, 1 ) ),
                  -yc,
                  0 ),
         fy * Z / L,
         hw + tf );
}

/// invalid input is rejected
void testInvalidInput()
{
  bool thrown = false;
  try {
    Beam beam( "B3D2",
               { 0., 0., 0., 0., 0., 10. },
               properties3D( rectanglePoints( -1, 1, -1, 1, 1, 1 ), { 0., 0., 1. }, 1. ),
               "LINEARELASTIC",
               linearElastic );
  }
  catch ( const std::invalid_argument& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, "orientation vector parallel to the axis must be rejected" );
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactoryNames,
    testCantilever2D,
    testCantilever3D,
    testUnsymmetricLProfile3D,
    testRigidBodyMotions,
    testBodyForce,
    testPlasticTangents,
    testPlasticLimits,
    testPlasticHingeCantilever,
    testInvalidInput,
  } );
  return 0;
}
