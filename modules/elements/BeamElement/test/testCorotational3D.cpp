#include "Marmot/CorotationalBeamElement3D.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Discretizations::ElementProperties;
using CR3D = Marmot::Discretizations::Elements::CorotationalBeamElement3D;

namespace {

  constexpr double pi = std::numbers::pi;

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

    Beam( Beam&& )            = default;
    Beam& operator=( Beam&& ) = delete;

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

    Eigen::VectorXd state( const std::string& name, int qp = 0 )
    {
      const auto view = element->getStateView( name, qp );
      return Eigen::Map< const Eigen::VectorXd >( view.stateLocation, view.stateSize );
    }
  };

  const double                E = 210000., nu = 0.3;
  const std::vector< double > linearElastic = { E, nu, 7.85e-9 };

  void check( double value, double expected, double tol, const std::string& what )
  {
    throwExceptionOnFailure( std::abs( value - expected ) <= tol,
                             Marmot::Utilities::MakeString() << std::setprecision( 12 ) << what << ": " << value
                                                             << " != " << expected << " (tol " << tol << ")" );
  }

  double checkTangent( Beam&                  beam,
                       const Eigen::VectorXd& U,
                       const Eigen::VectorXd& UOld,
                       double                 h,
                       double                 relTol,
                       const std::string&     what )
  {
    const auto [P, K]  = beam.kernels( U, UOld, false );
    const auto   numK  = beam.numericalTangent( U, UOld, h );
    const double scale = numK.cwiseAbs().maxCoeff();
    const double err   = ( K - numK ).cwiseAbs().maxCoeff() / scale;
    std::cout << what << ": tangent max rel. error " << err << "\n";
    throwExceptionOnFailure( err < relTol,
                             Marmot::Utilities::MakeString()
                               << what << ": tangent does not match the numerical one\nK=\n"
                               << K << "\nnumK=\n"
                               << numK );
    return err;
  }

  /// 3D section properties [vx, vy, vz, J, n, (y, z, A)...] of a rectangle h (along y) x b (along z), ny x nz cells
  /// with 2 x 2 Gauss points each (exact A, Iy, Iz); J given
  std::vector< double > rectangle3D( double h, double b, int ny, int nz, double J, Eigen::Vector3d v = { 0, 1, 0 } )
  {
    const double          g = 1 / std::sqrt( 3. ), dy = h / ny, dz = b / nz;
    std::vector< double > p = { v( 0 ), v( 1 ), v( 2 ), J, double( 4 * ny * nz ) };
    for ( int i = 0; i < ny; i++ )
      for ( int k = 0; k < nz; k++ )
        for ( double a : { -g, g } )
          for ( double c : { -g, g } )
            p.insert( p.end(),
                      { -h / 2 + ( i + 0.5 + 0.5 * a ) * dy, -b / 2 + ( k + 0.5 + 0.5 * c ) * dz, 0.25 * dy * dz } );
    return p;
  }

  /// the 2D section [n, (y, A)...] of the same rectangle (for B2D2CR, bending about z)
  std::vector< double > rectangle2D( double h, double b, int ny )
  {
    const double          g = 1 / std::sqrt( 3. ), dy = h / ny;
    std::vector< double > p = { double( 2 * ny ) };
    for ( int i = 0; i < ny; i++ )
      for ( double a : { -g, g } )
        p.insert( p.end(), { -h / 2 + ( i + 0.5 + 0.5 * a ) * dy, 0.5 * b * dy } );
    return p;
  }

  Eigen::Matrix3d rotation( const Eigen::Vector3d& w )
  {
    return CR3D::expSO3< double >( w );
  }
  Eigen::Vector3d rotationVector( const Eigen::Matrix3d& R )
  {
    return CR3D::logSO3< double >( R );
  }

  /// the global dofs of a rigid motion (rotation Q about X1, translation t) superposed on local dofs Ulocal (given in
  /// the reference configuration: displacements and rotation vectors relative to the reference)
  Eigen::VectorXd rigidMotion( const Eigen::Vector3d& X1,
                               const Eigen::Vector3d& X2,
                               const Eigen::VectorXd& Ulocal,
                               const Eigen::Matrix3d& Q,
                               const Eigen::Vector3d& t )
  {
    Eigen::VectorXd U( 12 );
    for ( int a = 0; a < 2; a++ ) {
      const Eigen::Vector3d X     = a == 0 ? X1 : X2;
      const Eigen::Vector3d x     = X1 + Q * ( X - X1 + Ulocal.segment< 3 >( 6 * a ) ) + t;
      U.segment< 3 >( 6 * a )     = x - X;
      U.segment< 3 >( 6 * a + 3 ) = rotationVector( Q * rotation( Ulocal.segment< 3 >( 6 * a + 3 ) ) );
    }
    return U;
  }

  /// A chain of beam elements through given nodes, clamped at the first node, Newton iterations under load control
  /// with dead nodal loads at the last node
  struct Chain {
    std::vector< Beam > beams;
    int                 nEl, nDof, nd;
    Eigen::VectorXd     U, UOld;
    Eigen::VectorXd     applied;

    Chain( const std::string&                    type,
           const std::vector< Eigen::VectorXd >& nodes,
           const std::vector< double >&          section,
           const std::string&                    material,
           const std::vector< double >&          matProps )
      : nEl( static_cast< int >( nodes.size() ) - 1 )
    {
      const int dim = static_cast< int >( nodes[0].size() );
      nd            = dim == 2 ? 3 : 6;
      nDof          = nd * ( nEl + 1 );
      for ( int e = 0; e < nEl; e++ ) {
        std::vector< double > X( nodes[e].data(), nodes[e].data() + dim );
        X.insert( X.end(), nodes[e + 1].data(), nodes[e + 1].data() + dim );
        beams.emplace_back( type, X, section, material, matProps );
      }
      U       = Eigen::VectorXd::Zero( nDof );
      UOld    = U;
      applied = Eigen::VectorXd::Zero( nd );
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > assemble( bool commit )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      for ( int e = 0; e < nEl; e++ ) {
        const auto [Pe, Ke] = beams[e].kernels( U.segment( nd * e, 2 * nd ), UOld.segment( nd * e, 2 * nd ), commit );
        P.segment( nd * e, 2 * nd ) += Pe;
        K.block( nd * e, nd * e, 2 * nd, 2 * nd ) += Ke;
      }
      return { P, K };
    }

    int solve( const Eigen::VectorXd& tipLoad, int nSteps, double tol = 1e-10 )
    {
      int                   total = 0;
      Eigen::VectorXd       F     = Eigen::VectorXd::Zero( nDof );
      const Eigen::VectorXd start = applied;
      for ( int step = 1; step <= nSteps; step++ ) {
        F.tail( nd )   = start + ( tipLoad - start ) * double( step ) / nSteps;
        bool converged = false;
        for ( int it = 0; it < 50; it++, total++ ) {
          const auto [P, K]       = assemble( false );
          const Eigen::VectorXd r = ( P - F ).tail( nDof - nd );
          if ( r.norm() < tol * std::max( 1.0, F.norm() ) ) {
            converged = true;
            break;
          }
          U.tail( nDof - nd ) += K.bottomRightCorner( nDof - nd, nDof - nd ).lu().solve( -r );
        }
        throwExceptionOnFailure( converged,
                                 Marmot::Utilities::MakeString() << "chain: Newton did not converge in step " << step );
        assemble( true );
        UOld    = U;
        applied = F.tail( nd );
      }
      return total;
    }

    Eigen::VectorXd tip() const { return U.tail( nd ); }
  };

  std::vector< Eigen::VectorXd > straightNodes( int nEl, double L, int dim )
  {
    std::vector< Eigen::VectorXd > n;
    for ( int e = 0; e <= nEl; e++ ) {
      Eigen::VectorXd x = Eigen::VectorXd::Zero( dim );
      x( 0 )            = e * L / nEl;
      n.push_back( x );
    }
    return n;
  }

  /// elastica of a cantilever under a dead tip force normal to the axis (as in testCorotational.cpp)
  std::array< double, 3 > elastica( double lambda )
  {
    const int n     = 20000;
    auto      shoot = [&]( double thetaTip, double& x, double& y, double& th ) {
      double t = 0, k = std::sqrt( 2 * lambda * std::sin( thetaTip ) );
      x = y          = 0;
      const double h = 1.0 / n;
      auto         f = [&]( double t_, double k_ ) {
        return std::array< double, 4 >{ k_, -lambda * std::cos( t_ ), std::cos( t_ ), std::sin( t_ ) };
      };
      for ( int i = 0; i < n; i++ ) {
        const auto a = f( t, k );
        const auto b = f( t + 0.5 * h * a[0], k + 0.5 * h * a[1] );
        const auto c = f( t + 0.5 * h * b[0], k + 0.5 * h * b[1] );
        const auto d = f( t + h * c[0], k + h * c[1] );
        t += h / 6 * ( a[0] + 2 * b[0] + 2 * c[0] + d[0] );
        k += h / 6 * ( a[1] + 2 * b[1] + 2 * c[1] + d[1] );
        x += h / 6 * ( a[2] + 2 * b[2] + 2 * c[2] + d[2] );
        y += h / 6 * ( a[3] + 2 * b[3] + 2 * c[3] + d[3] );
      }
      th = t;
      return k;
    };
    double x, y, th, lo = 1e-9, hi = pi / 2 - 1e-12;
    for ( int it = 0; it < 60; it++ ) {
      const double mid                         = 0.5 * ( lo + hi );
      ( shoot( mid, x, y, th ) < 0 ? lo : hi ) = mid;
    }
    shoot( 0.5 * ( lo + hi ), x, y, th );
    return { y, 1 - x, th };
  }

} // namespace

/// B3D2CR is registered and has the B3D2 interface; exp/log are inverse
void testFactoryAndRotations()
{
  std::unique_ptr< Marmot::Discretizations::Elements::Element > el(
    Marmot::Registry::ElementFactory::createElement( "B3D2CR", 1 ) );
  throwExceptionOnFailure( el->getNSpatialDimensions() == 3 && el->getNDofPerElement() == 12 &&
                             el->getElementShape() == "bar2" && el->getNNodes() == 2 &&
                             el->getNodeFields()[1] == std::vector< std::string >{ "displacement", "rotation" },
                           "B3D2CR element description" );
  for ( const Eigen::Vector3d& w : { Eigen::Vector3d( 0, 0, 0 ),
                                     Eigen::Vector3d( 1e-7, -2e-7, 3e-8 ),
                                     Eigen::Vector3d( 0.3, -1.1, 0.7 ),
                                     Eigen::Vector3d( Eigen::Vector3d( 2.0, 1.5, -1.0 ).normalized() * 3.0 ) } ) {
    const Eigen::Matrix3d R = rotation( w );
    throwExceptionOnFailure( ( R.transpose() * R - Eigen::Matrix3d::Identity() ).cwiseAbs().maxCoeff() < 1e-14 &&
                               ( rotationVector( R ) - w ).cwiseAbs().maxCoeff() < 1e-12,
                             Marmot::Utilities::MakeString() << "exp/log of " << w.transpose() );
  }
}

/// rigid rotations up to 170 degrees about inclined axes (and a translation): zero internal forces; superposed on a
/// deformed state: the nodal forces rotate with the beam and the section forces are unchanged
void testRigidBodyRotations()
{
  const Eigen::Vector3d X1( 10., 20., -5. ), X2( 70., 100., 15. );
  const auto            section = rectangle3D( 10, 6, 1, 1, 400., { 0, 0, 1 } );
  Eigen::VectorXd       Udef( 12 );
  Udef << 0., 0., 0., 0.01, -0.02, 0.015, 0.05, -0.3, 0.2, -0.02, 0.03, -0.01;

  Beam                  ref( "B3D2CR",
                             { X1( 0 ), X1( 1 ), X1( 2 ), X2( 0 ), X2( 1 ), X2( 2 ) },
            section,
            "LINEARELASTIC",
            linearElastic );
  const Eigen::VectorXd zero  = Eigen::VectorXd::Zero( 12 );
  const auto [P0, K0]         = ref.kernels( Udef, zero, true );
  const Eigen::VectorXd s0    = ref.state( "section forces", 1 );
  const double          scale = P0.cwiseAbs().maxCoeff();

  for ( const auto& [deg, axis] : std::vector< std::pair< double, Eigen::Vector3d > >{ { 90., { 0, 0, 1 } },
                                                                                       { 170., { 1, 1, 0 } },
                                                                                       { 60., { -1, 2, 3 } },
                                                                                       { 135., { 0.3, -1, 0.2 } } } ) {
    const Eigen::Matrix3d Q = rotation( deg * pi / 180 * axis.normalized() );
    Beam                  beam( "B3D2CR",
                                { X1( 0 ), X1( 1 ), X1( 2 ), X2( 0 ), X2( 1 ), X2( 2 ) },
               section,
               "LINEARELASTIC",
               linearElastic );
    const auto            Urigid = rigidMotion( X1, X2, zero, Q, { 3., -7., 2. } );
    const auto [P, K]            = beam.kernels( Urigid, zero, false );
    std::cout << "rigid rotation " << deg << " deg: |P| = " << P.cwiseAbs().maxCoeff() << "\n";
    throwExceptionOnFailure( P.cwiseAbs().maxCoeff() < 1e-11 * scale,
                             Marmot::Utilities::MakeString()
                               << "rigid rotation " << deg << ": internal forces " << P.transpose() );

    const auto Ud       = rigidMotion( X1, X2, Udef, Q, { 3., -7., 2. } );
    const auto [Pd, Kd] = beam.kernels( Ud, zero, false );
    double dF           = 0;
    for ( int a = 0; a < 2; a++ )
      dF = std::max( dF, ( Pd.segment< 3 >( 6 * a ) - Q * P0.segment< 3 >( 6 * a ) ).cwiseAbs().maxCoeff() );
    throwExceptionOnFailure( dF < 1e-9 * scale,
                             Marmot::Utilities::MakeString() << "objectivity " << deg << ": nodal forces" );
    beam.kernels( Ud, zero, true );
    throwExceptionOnFailure( ( beam.state( "section forces", 1 ) - s0 ).cwiseAbs().maxCoeff() <
                               1e-9 * s0.cwiseAbs().maxCoeff(),
                             "objectivity: section forces" );
    checkTangent( beam, Ud, zero, 1e-7, 1e-6, Marmot::Utilities::MakeString() << "elastic rotated " << deg );
  }
}

/// small deformations: K(U = 0) of B3D2CR equals B3D2; forces agree to first order in the rotations
void testSmallDeformationLimit()
{
  const std::vector< double > X       = { 0., 0., 0., 60., 80., 30. };
  const auto                  section = rectangle3D( 10, 6, 2, 1, 400., { 0, 0, 1 } );
  Beam                        lin( "B3D2", X, section, "LINEARELASTIC", linearElastic );
  Beam                        cr( "B3D2CR", X, section, "LINEARELASTIC", linearElastic );
  const Eigen::VectorXd       zero = Eigen::VectorXd::Zero( 12 );
  const auto                  Klin = lin.kernels( zero, zero, false ).second;
  const auto                  Kcr  = cr.kernels( zero, zero, false ).second;
  const double                dK   = ( Klin - Kcr ).cwiseAbs().maxCoeff() / Klin.cwiseAbs().maxCoeff();
  std::cout << "K(U=0): max rel. difference B3D2CR - B3D2 = " << dK << "\n";
  throwExceptionOnFailure( dK < 1e-12, "small deformation limit: tangent at U = 0" );

  Eigen::VectorXd dir( 12 );
  dir << 0.1, -0.2, 0.05, 0.003, -0.002, 0.001, 0.3, 0.1, -0.1, -0.002, 0.001, 0.004;
  for ( double eps : { 1e-2, 1e-3, 1e-4 } ) {
    const Eigen::VectorXd U    = eps * dir;
    const auto            Plin = lin.kernels( U, zero, false ).first;
    const auto            Pcr  = cr.kernels( U, zero, false ).first;
    const double          d    = ( Plin - Pcr ).norm() / Plin.norm();
    std::cout << "small U (scale " << eps << "): |P_CR - P_lin| / |P_lin| = " << d << "\n";
    throwExceptionOnFailure( d < 10 * eps, "small deformation limit: forces" );
  }
}

/// numerical tangents: elastic and von Mises (plastic) in large rotation states
void testTangents()
{
  const std::vector< double > vm = { 210000., 0.3, 500., 2000., 0., 1., 7.85e-9 };
  const Eigen::Vector3d       X1( 0., 0., 0. ), X2( 60., 80., -20. );
  const auto                  section = rectangle3D( 20, 10, 3, 2, 2800., { 0, 0, 1 } );
  {
    Beam            beam( "B3D2CR", { 0., 0., 0., 60., 80., -20. }, section, "LINEARELASTIC", linearElastic );
    Eigen::VectorXd Ul( 12 );
    Ul << 0., 0., 0., 0.05, -0.03, 0.04, 0.08, -2.0, 1.0, 0.02, 0.06, -0.05;
    checkTangent( beam,
                  rigidMotion( X1, X2, Ul, rotation( { 0.4, -0.9, 0.6 } ), { 5., 1., -2. } ),
                  Eigen::VectorXd::Zero( 12 ),
                  1e-7,
                  1e-6,
                  "elastic, 3D rotation" );
  }
  {
    Beam            beam( "B3D2CR", { 0., 0., 0., 60., 80., -20. }, section, "VONMISES", vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 12 ), U;
    auto            path = [&]( double f ) {
      Eigen::VectorXd Ul( 12 );
      Ul << 0., 0., 0., 0.02 * f, -0.03 * f, 0.01 * f, 0.25 * f, -0.6 * f, 0.5 * f, -0.01 * f, 0.04 * f, 0.03 * f;
      return rigidMotion( X1, X2, Ul, rotation( f * Eigen::Vector3d( 0.5, 0.8, -0.4 ) ), { 2. * f, -3 * f, f } );
    };
    for ( int i = 1; i <= 10; i++ ) {
      U = path( i / 10. );
      beam.kernels( U, UOld, true );
      UOld = U;
    }
    double kappaMax = 0;
    for ( int f = 0; f < 24; f++ )
      kappaMax = std::max( kappaMax,
                           beam.state( Marmot::Utilities::MakeString() << "fiber " << f << " kappa", 0 )( 0 ) );
    std::cout << "von Mises path: max fiber kappa at QP 0 = " << kappaMax << "\n";
    throwExceptionOnFailure( kappaMax > 1e-4, "von Mises: plastic state" );
    checkTangent( beam, path( 1.1 ), UOld, 1e-7, 1e-5, "von Mises, plastic step" );
  }
}

/// in-plane problems equal B2D2CR: roll-up under an end moment (to 3/4 of a circle; the rotation vector
/// parametrization is singular at 2 pi) and the elastica under a dead tip force
void testInPlaneEqualsB2D2CR()
{
  const double L = 100, h = 2, b = 1;
  const double Iz = b * h * h * h / 12, EI = E * Iz;
  const auto   sec3 = rectangle3D( h, b, 2, 1, 1., { 0, 1, 0 } );
  const auto   sec2 = rectangle2D( h, b, 2 );
  {
    const int nEl = 10;
    Chain     c3( "B3D2CR", straightNodes( nEl, L, 3 ), sec3, "LINEARELASTIC", linearElastic );
    Chain     c2( "B2D2CR", straightNodes( nEl, L, 2 ), sec2, "LINEARELASTIC", linearElastic );
    for ( int part = 1; part <= 3; part++ ) {
      const double    M = part * 0.25 * 2 * pi * EI / L, R = EI / M;
      Eigen::VectorXd F3 = Eigen::VectorXd::Zero( 6 ), F2 = Eigen::VectorXd::Zero( 3 );
      F3( 5 ) = M, F2( 2 ) = M;
      const int its = c3.solve( F3, 5 );
      c2.solve( F2, 5 );
      const auto   t3 = c3.tip(), t2 = c2.tip();
      const double x = L + t3( 0 ), y = t3( 1 );
      const double err = std::hypot( x - R * std::sin( L / R ), y - R * ( 1 - std::cos( L / R ) ) ) / L;
      std::cout << std::setprecision( 9 ) << "roll-up M = " << part * 0.25 << " x 2 pi EI/L: tip (" << x << ", " << y
                << ", " << t3( 2 ) << "), rotation " << t3( 5 ) << " (exact " << L / R << "), vs B2D2CR "
                << ( t3.head< 2 >() - t2.head< 2 >() ).norm() / L << ", position error / L " << err << ", " << its
                << " iterations\n";
      check( t3( 5 ), L / R, 1e-8 * L / R, "roll-up: tip rotation" );
      check( ( t3.head< 2 >() - t2.head< 2 >() ).norm() / L, 0, 1e-9, "roll-up: B3D2CR vs B2D2CR" );
      check( t3( 2 ), 0, 1e-10, "roll-up: out of plane" );
    }
  }
  {
    const auto   s3  = rectangle3D( 1, 1, 1, 1, 1., { 0, 1, 0 } );
    const auto   s2  = rectangle2D( 1, 1, 1 );
    const double EI1 = E / 12;
    Chain        c3( "B3D2CR", straightNodes( 20, L, 3 ), s3, "LINEARELASTIC", linearElastic );
    Chain        c2( "B2D2CR", straightNodes( 20, L, 2 ), s2, "LINEARELASTIC", linearElastic );
    for ( double lambda : { 1., 5., 10. } ) {
      const double    P  = lambda * EI1 / ( L * L );
      const auto      rf = elastica( lambda );
      Eigen::VectorXd F3 = Eigen::VectorXd::Zero( 6 ), F2 = Eigen::VectorXd::Zero( 3 );
      F3( 1 ) = P, F2( 1 ) = P;
      c3.solve( F3, 10, 1e-8 );
      c2.solve( F2, 10, 1e-8 );
      const auto t3 = c3.tip(), t2 = c2.tip();
      std::cout << std::setprecision( 6 ) << "elastica lambda " << lambda << ": v/L " << t3( 1 ) / L << " (" << rf[0]
                << "), (L-x)/L " << -t3( 0 ) / L << " (" << rf[1] << "), theta " << t3( 5 ) << " (" << rf[2]
                << "), |B3D2CR - B2D2CR| / L " << ( t3.head< 2 >() - t2.head< 2 >() ).norm() / L << "\n";
      check( t3( 1 ) / L, rf[0], 1e-3, "elastica v/L" );
      check( -t3( 0 ) / L, rf[1], 1e-3, "elastica (L-x)/L" );
      check( t3( 5 ), rf[2], 1e-3, "elastica theta" );
      check( ( t3.head< 2 >() - t2.head< 2 >() ).norm() / L, 0, 1e-8, "elastica: B3D2CR vs B2D2CR" );
    }
  }
}

/// 45 degree bend (Bathe and Bolourchi 1979): cantilever arc of radius 100 in the x-y plane, unit square section,
/// E = 1e7, nu = 0, out-of-plane dead tip load P = 300, 450, 600; 8 elements. Reference tip displacements (Simo and
/// Vu-Quoc 1986, 8 elements; Crisfield 1990): P = 600: (u, v, w) = (-23.5, -13.4, 53.4) (Bathe), (-23.7, -13.6, 53.6)
void testBend45()
{
  const double                   R = 100;
  std::vector< Eigen::VectorXd > nodes;
  for ( int i = 0; i <= 8; i++ ) {
    const double    a = pi / 4 * i / 8;
    Eigen::VectorXd x( 3 );
    x << R * std::sin( a ), R * ( 1 - std::cos( a ) ), 0.;
    nodes.push_back( x );
  }
  // unit square, J = 0.1406 (Saint-Venant torsion constant of a square), the section points carry it via the scaling
  const auto                  section = rectangle3D( 1, 1, 4, 4, 0.1406, { 0, 0, 1 } );
  const std::vector< double > mat     = { 1e7, 0.0, 1. };
  Chain                       chain( "B3D2CR", nodes, section, "LINEARELASTIC", mat );
  int                         its = 0;
  std::array< double, 3 >     ref600{ -23.5, -13.4, 53.4 };
  for ( double P : { 300., 450., 600. } ) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero( 6 );
    F( 2 )            = P;
    its += chain.solve( F, 5, 1e-9 );
    const auto t = chain.tip();
    // tip displacements in the frame of Bathe: u along the tangent direction... reported as global x, y, z changes
    std::cout << std::setprecision( 6 ) << "45 deg bend P = " << P << ": tip displacement (" << t( 0 ) << ", " << t( 1 )
              << ", " << t( 2 ) << ")\n";
    if ( P == 600 ) {
      check( t( 0 ), ref600[0], 0.03 * std::abs( ref600[0] ) + 0.5, "45 deg bend: u" );
      check( t( 1 ), ref600[1], 0.03 * std::abs( ref600[1] ) + 0.5, "45 deg bend: v" );
      check( t( 2 ), ref600[2], 0.02 * std::abs( ref600[2] ), "45 deg bend: w" );
    }
  }
  std::cout << "  Newton iterations: " << its << "\n";
}

/// von Mises fiber section: pure bending about y (moment-curvature path, loading and unloading) is independent of a
/// superposed large rigid 3D rotation and equals B3D2 when the chord does not rotate
void testPlasticBending()
{
  const double                L = 100, fy = 500, hgt = 20, b = 10;
  const auto                  section = rectangle3D( hgt, b, 1, 6, 2800., { 0, 1, 0 } );
  const std::vector< double > vm      = { E, nu, fy, 0., 0., 1., 7.85e-9 };
  const double                Mp      = fy * hgt * b * b / 4; // bending about y: lever along z (width b)
  const double                kappaY  = 2 * fy / ( E * b );
  Beam                        lin( "B3D2", { 0., 0., 0., L, 0., 0. }, section, "VONMISES", vm );
  Beam                        cr( "B3D2CR", { 0., 0., 0., L, 0., 0. }, section, "VONMISES", vm );
  Beam                        crRot( "B3D2CR", { 0., 0., 0., L, 0., 0. }, section, "VONMISES", vm );
  Eigen::VectorXd             UOld = Eigen::VectorXd::Zero( 12 ), URotOld = UOld;
  double                      maxDiff = 0, maxDiffRot = 0, Mmax = 0;
  for ( int i = 1; i <= 40; i++ ) {
    const double    k     = i <= 20 ? 6 * kappaY * i / 20. : 6 * kappaY - 9 * kappaY * ( i - 20 ) / 20.;
    const double    Theta = 0.5 * k * L;
    Eigen::VectorXd U     = Eigen::VectorXd::Zero( 12 );
    U( 4 ) = -Theta, U( 10 ) = Theta; // rotations about y
    lin.kernels( U, UOld, true );
    cr.kernels( U, UOld, true );
    const auto URot = rigidMotion( { 0., 0., 0. },
                                   { L, 0., 0. },
                                   U,
                                   rotation( ( 1.0 + 0.02 * i ) * Eigen::Vector3d( 0.6, -0.3, 0.74 ) ),
                                   { 1., 2., 3. } );
    crRot.kernels( URot, URotOld, true );
    UOld = U, URotOld = URot;
    const double Ml = lin.state( "bending moments", 1 )( 0 ), Mc = cr.state( "bending moments", 1 )( 0 ),
                 Mr = crRot.state( "bending moments", 1 )( 0 );
    maxDiff         = std::max( maxDiff, std::abs( Ml - Mc ) / Mp );
    maxDiffRot      = std::max( maxDiffRot, std::abs( Mr - Mc ) / Mp );
    Mmax            = std::max( Mmax, std::abs( Mc ) / Mp );
  }
  std::cout << "plastic bending: max |M_CR - M_B3D2| / M_p = " << maxDiff
            << ", rigidly rotating CR vs CR: " << maxDiffRot << ", max M / M_p = " << Mmax << "\n";
  throwExceptionOnFailure( maxDiff < 1e-10 && maxDiffRot < 1e-8 && Mmax > 0.9, "plastic bending: CR = B3D2" );
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactoryAndRotations,
    testRigidBodyRotations,
    testSmallDeformationLimit,
    testTangents,
    testInPlaneEqualsB2D2CR,
    testBend45,
    testPlasticBending,
  } );
  return 0;
}
