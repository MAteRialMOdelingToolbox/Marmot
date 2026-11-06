#include "Marmot/CorotationalBeamElement.h"
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

    /// residual and tangent for the total dofs U from the committed state at UOld; commits if requested
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

  /// section points (y, A) of a rectangle of height h, width b: n layers with 2 Gauss points each (exact A, I)
  std::vector< double > rectangle2D( double h, double b, int n )
  {
    const double          g = 1 / std::sqrt( 3. ), dy = h / n;
    std::vector< double > p = { double( 2 * n ) };
    for ( int i = 0; i < n; i++ )
      for ( double a : { -g, g } )
        p.insert( p.end(), { -h / 2 + ( i + 0.5 + 0.5 * a ) * dy, 0.5 * b * dy } );
    return p;
  }

  /// the global dofs of a rigid body motion (translation t, rotation phi about the first node) of a beam X1 -> X2
  /// superposed on local dofs (relative to the reference) Ulocal, given in the reference configuration
  Eigen::VectorXd rigidMotion( const Eigen::Vector2d& X1,
                               const Eigen::Vector2d& X2,
                               const Eigen::VectorXd& Ulocal,
                               double                 phi,
                               const Eigen::Vector2d& t )
  {
    const Eigen::Rotation2Dd R( phi );
    Eigen::VectorXd          U( 6 );
    for ( int a = 0; a < 2; a++ ) {
      const Eigen::Vector2d X = a == 0 ? X1 : X2;
      const Eigen::Vector2d x = X1 + R * ( X - X1 + Ulocal.segment< 2 >( 3 * a ) ) + t;
      U.segment< 2 >( 3 * a ) = x - X;
      U( 3 * a + 2 )          = Ulocal( 3 * a + 2 ) + phi;
    }
    return U;
  }

  /// A straight chain of beam elements from (0, 0) along x, clamped at the first node, solved by Newton iterations on
  /// the assembled system under load control with dead tip loads
  struct Chain {
    std::vector< Beam > beams;
    int                 nEl, nDof;
    Eigen::VectorXd     U, UOld;
    Eigen::Vector3d     applied = Eigen::Vector3d::Zero(); ///< the tip load of the last converged state

    Chain( const std::string&           type,
           int                          nElements,
           double                       L,
           const std::vector< double >& section,
           const std::string&           material,
           const std::vector< double >& matProps )
      : nEl( nElements ), nDof( 3 * ( nElements + 1 ) )
    {
      for ( int e = 0; e < nEl; e++ )
        beams.emplace_back( type,
                            std::vector< double >{ e * L / nEl, 0., ( e + 1 ) * L / nEl, 0. },
                            section,
                            material,
                            matProps );
      U    = Eigen::VectorXd::Zero( nDof );
      UOld = U;
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > assemble( bool commit )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      for ( int e = 0; e < nEl; e++ ) {
        const auto [Pe, Ke] = beams[e].kernels( U.segment( 3 * e, 6 ), UOld.segment( 3 * e, 6 ), commit );
        P.segment( 3 * e, 6 ) += Pe;
        K.block( 3 * e, 3 * e, 6, 6 ) += Ke;
      }
      return { P, K };
    }

    /// increase the tip load [Fx, Fy, M] (dead) from the current one to tipLoad in nSteps equal increments; returns
    /// the number of Newton iterations
    int solve( const Eigen::Vector3d& tipLoad, int nSteps, double tol = 1e-10 )
    {
      int                   total = 0;
      Eigen::VectorXd       F     = Eigen::VectorXd::Zero( nDof );
      const Eigen::Vector3d start = applied;
      for ( int step = 1; step <= nSteps; step++ ) {
        F.tail< 3 >()  = start + ( tipLoad - start ) * double( step ) / nSteps;
        bool converged = false;
        for ( int it = 0; it < 50; it++, total++ ) {
          const auto [P, K]       = assemble( false );
          const Eigen::VectorXd r = ( P - F ).tail( nDof - 3 );
          if ( r.norm() < tol * std::max( 1.0, F.norm() ) ) {
            converged = true;
            break;
          }
          U.tail( nDof - 3 ) += K.bottomRightCorner( nDof - 3, nDof - 3 ).lu().solve( -r );
        }
        throwExceptionOnFailure( converged,
                                 Marmot::Utilities::MakeString() << "chain: Newton did not converge in step " << step );
        assemble( true );
        UOld    = U;
        applied = F.tail< 3 >();
      }
      return total;
    }

    Eigen::Vector3d tip() const { return U.tail< 3 >(); }
  };

  /// the elastica of a cantilever under a dead tip force P normal to the undeformed axis: tip deflection v/L, axial
  /// shortening (L - x)/L and tip rotation for lambda = P L^2 / EI. theta'' = -lambda cos(theta) on s in [0, 1],
  /// theta(0) = 0, theta'(1) = 0; first integral theta'(0)^2 = 2 lambda sin(theta_tip): bisection on theta_tip (RK4)
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
      return k; // theta'(1): < 0 if theta_tip is too small, > 0 if too large
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

/// B2D2CR is registered and has the B2D2 interface
void testFactory()
{
  std::unique_ptr< Marmot::Discretizations::Elements::Element > el(
    Marmot::Registry::ElementFactory::createElement( "B2D2CR", 1 ) );
  throwExceptionOnFailure( el->getNSpatialDimensions() == 2 && el->getNDofPerElement() == 6 &&
                             el->getElementShape() == "bar2" && el->getNNodes() == 2 &&
                             el->getNodeFields()[1] == std::vector< std::string >{ "displacement", "rotation" },
                           "B2D2CR element description" );
}

/// rigid body rotations of 90, 180, 270 and 540 degrees (and a translation) of an undeformed inclined beam: zero
/// internal forces; superposed on a deformed state: the forces rotate with the beam (objectivity) and the section
/// forces are unchanged
void testRigidBodyRotations()
{
  const Eigen::Vector2d X1( 10., 20. ), X2( 70., 100. ); // L = 100
  const auto            section = rectangle2D( 10, 5, 2 );
  Eigen::VectorXd       Udef( 6 );
  Udef << 0., 0., 0.01, 0.05, -0.3, -0.02; // a deformed state (local), relative to the reference

  Beam ref( "B2D2CR", { X1( 0 ), X1( 1 ), X2( 0 ), X2( 1 ) }, section, "LINEARELASTIC", linearElastic );
  const auto [P0, K0]         = ref.kernels( Udef, Eigen::VectorXd::Zero( 6 ), true );
  const Eigen::VectorXd s0    = ref.state( "section forces", 1 );
  const double          scale = P0.cwiseAbs().maxCoeff();

  for ( double deg : { 90., 180., 270., 540., -135. } ) {
    const double phi = deg * pi / 180;
    Beam         beam( "B2D2CR", { X1( 0 ), X1( 1 ), X2( 0 ), X2( 1 ) }, section, "LINEARELASTIC", linearElastic );
    const auto   Urigid = rigidMotion( X1, X2, Eigen::VectorXd::Zero( 6 ), phi, { 3., -7. } );
    const auto [P, K]   = beam.kernels( Urigid, Eigen::VectorXd::Zero( 6 ), false );
    std::cout << "rigid rotation " << deg << " deg: |P| = " << P.cwiseAbs().maxCoeff() << "\n";
    throwExceptionOnFailure( P.cwiseAbs().maxCoeff() < 1e-12 * scale,
                             Marmot::Utilities::MakeString()
                               << "rigid rotation " << deg << ": internal forces " << P.transpose() );

    // deformed + rotated: the global forces rotate, the section forces are those of the unrotated state
    const auto Ud       = rigidMotion( X1, X2, Udef, phi, { 3., -7. } );
    const auto [Pd, Kd] = beam.kernels( Ud, Eigen::VectorXd::Zero( 6 ), false );
    const Eigen::Rotation2Dd R( phi );
    Eigen::VectorXd          Pexpected( 6 );
    for ( int a = 0; a < 2; a++ ) {
      Pexpected.segment< 2 >( 3 * a ) = R * P0.segment< 2 >( 3 * a );
      Pexpected( 3 * a + 2 )          = P0( 3 * a + 2 );
    }
    throwExceptionOnFailure( ( Pd - Pexpected ).cwiseAbs().maxCoeff() < 1e-9 * scale,
                             Marmot::Utilities::MakeString()
                               << "objectivity " << deg << ": " << Pd.transpose() << " vs " << Pexpected.transpose() );
    beam.kernels( Ud, Eigen::VectorXd::Zero( 6 ), true );
    throwExceptionOnFailure( ( beam.state( "section forces", 1 ) - s0 ).cwiseAbs().maxCoeff() <
                               1e-9 * s0.cwiseAbs().maxCoeff(),
                             "objectivity: section forces" );
    // tangent in the rotated deformed state
    checkTangent( beam,
                  Ud,
                  Eigen::VectorXd::Zero( 6 ),
                  1e-7,
                  1e-6,
                  Marmot::Utilities::MakeString() << "elastic rotated " << deg );
  }
}

/// small deformations: B2D2CR = B2D2 (tangent at U = 0 identical; forces for a small U to second order)
void testSmallDeformationLimit()
{
  const std::vector< double > X       = { 0., 0., 60., 80. };
  const auto                  section = rectangle2D( 10, 5, 2 );
  Beam                        lin( "B2D2", X, section, "LINEARELASTIC", linearElastic );
  Beam                        cr( "B2D2CR", X, section, "LINEARELASTIC", linearElastic );
  const Eigen::VectorXd       zero = Eigen::VectorXd::Zero( 6 );
  const auto                  Klin = lin.kernels( zero, zero, false ).second;
  const auto                  Kcr  = cr.kernels( zero, zero, false ).second;
  const double                dK   = ( Klin - Kcr ).cwiseAbs().maxCoeff() / Klin.cwiseAbs().maxCoeff();
  std::cout << "K(U=0): max rel. difference B2D2CR - B2D2 = " << dK << "\n";
  throwExceptionOnFailure( dK < 1e-13, "small deformation limit: tangent at U = 0" );

  Eigen::VectorXd dir( 6 );
  dir << 0.1, -0.2, 0.003, 0.3, 0.1, -0.002;
  for ( double eps : { 1e-2, 1e-3, 1e-4 } ) {
    const Eigen::VectorXd U    = eps * dir;
    const auto            Plin = lin.kernels( U, zero, false ).first;
    const auto            Pcr  = cr.kernels( U, zero, false ).first;
    const double          d    = ( Plin - Pcr ).norm() / Plin.norm();
    std::cout << "small U (scale " << eps << "): |P_CR - P_lin| / |P_lin| = " << d << "\n";
    throwExceptionOnFailure( d < 10 * eps, "small deformation limit: forces" ); // first order in the rotations
  }
}

/// numerical tangents: elastic and von Mises in large rotation states
void testTangents()
{
  //                                 E        nu   fy    H     dfy  delta rho
  const std::vector< double > vm = { 210000., 0.3, 500., 2000., 0., 1., 7.85e-9 };
  const Eigen::Vector2d       X1( 0., 0. ), X2( 60., 80. );
  const auto                  section = rectangle2D( 20, 10, 3 );

  // elastic, large rotation with bending and stretching
  {
    Beam            beam( "B2D2CR", { 0., 0., 60., 80. }, section, "LINEARELASTIC", linearElastic );
    Eigen::VectorXd Ul( 6 );
    Ul << 0., 0., 0.05, 0.08, -2.0, 0.02;
    checkTangent( beam,
                  rigidMotion( X1, X2, Ul, 1.2, { 5., 1. } ),
                  Eigen::VectorXd::Zero( 6 ),
                  1e-7,
                  1e-6,
                  "elastic, 69 deg" );
  }
  // plastic: a path of bending + tension while the beam rotates by 70 degrees, then the tangent in a plastic step
  {
    Beam            beam( "B2D2CR", { 0., 0., 60., 80. }, section, "VONMISES", vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 6 ), U;
    for ( int i = 1; i <= 10; i++ ) {
      Eigen::VectorXd Ul( 6 );
      const double    f = i / 10.;
      // axis (0.6, 0.8), normal (-0.8, 0.6): elongation 0.4 (strain 4e-3), deflection 0.5 of the second node
      Ul << 0., 0., -0.02 * f, ( 0.24 - 0.4 ) * f, ( 0.32 + 0.3 ) * f, 0.03 * f;
      U = rigidMotion( X1, X2, Ul, 70. * pi / 180 * f, { 2. * f, -3 * f } );
      beam.kernels( U, UOld, true );
      UOld = U;
    }
    std::cout << "von Mises path: kappa of the section points at QP 0:";
    double kappaMax = 0;
    for ( int f = 0; f < 6; f++ ) {
      const double k = beam.state( Marmot::Utilities::MakeString() << "fiber " << f << " kappa", 0 )( 0 );
      std::cout << " " << k;
      kappaMax = std::max( kappaMax, k );
    }
    std::cout << "\n";
    throwExceptionOnFailure( kappaMax > 1e-4, "von Mises: plastic state" );
    Eigen::VectorXd Ul( 6 );
    Ul << 0., 0., -0.022, 0.252 - 0.4, 0.336 + 0.3, 0.033;
    U = rigidMotion( X1, X2, Ul, 75. * pi / 180, { 2.2, -3.3 } );
    checkTangent( beam, U, UOld, 1e-7, 1e-5, "von Mises, plastic step at 70 -> 75 deg" );
  }
}

/// end moment M = 2 pi EI / L rolls the cantilever up into a full circle: the tip at the root, the tip rotation 2 pi;
/// intermediate loads vs the circle x = R sin(L/R), y = R (1 - cos(L/R)), R = EI / M
void testRollUp()
{
  const double L = 100, h = 2, b = 1;
  const auto   section = rectangle2D( h, b, 2 );
  const double EI      = E * b * h * h * h / 12;
  for ( int nEl : { 10, 20, 40 } ) {
    Chain        chain( "B2D2CR", nEl, L, section, "LINEARELASTIC", linearElastic );
    const double Mfull   = 2 * pi * EI / L;
    double       errPrev = 0;
    int          its     = 0;
    for ( int part = 1; part <= 4; part++ ) {
      its += chain.solve( { 0, 0, part * 0.25 * Mfull }, 5 );
      const double M = part * 0.25 * Mfull, R = EI / M;
      const auto   tip = chain.tip();
      const double x = L + tip( 0 ), y = tip( 1 );
      const double err = std::hypot( x - R * std::sin( L / R ), y - R * ( 1 - std::cos( L / R ) ) ) / L;
      std::cout << std::setprecision( 8 ) << "roll-up " << nEl << " elements, M = " << part * 0.25
                << " x 2 pi EI/L: tip (" << x << ", " << y << "), rotation " << tip( 2 ) << " (exact " << L / R
                << "), position error / L = " << err << "\n";
      check( tip( 2 ), L / R, 1e-9 * L / R, "roll-up: tip rotation" );
      errPrev = err;
    }
    std::cout << "  Newton iterations: " << its << "\n";
    // full circle: the tip at the root; the chord of an element of a circle of radius R: error O(h^2)
    throwExceptionOnFailure( errPrev < 0.5 * ( pi / nEl ) * ( pi / nEl ), "roll-up: tip position (full circle)" );
  }
}

/// elastica: cantilever with a dead tip force normal to the axis, lambda = P L^2 / EI = 1, 2, 5, 10 vs the shooting
/// solution (Mattiasson's values)
void testElastica()
{
  const double L = 100, h = 1, b = 1;
  const auto   section = rectangle2D( h, b, 1 );
  const double EI      = E * b * h * h * h / 12;
  Chain        chain( "B2D2CR", 40, L, section, "LINEARELASTIC", linearElastic );
  int          its = 0;
  for ( double lambda : { 1., 2., 5., 10. } ) {
    const double P   = lambda * EI / ( L * L );
    const auto   ref = elastica( lambda );
    its += chain.solve( { 0, P, 0 }, 10, 1e-8 );
    const auto tip = chain.tip();
    std::cout << std::setprecision( 6 ) << "elastica lambda " << lambda << ": v/L " << tip( 1 ) / L << " (" << ref[0]
              << "), (L-x)/L " << -tip( 0 ) / L << " (" << ref[1] << "), theta " << tip( 2 ) << " (" << ref[2] << ")\n";
    check( tip( 1 ) / L, ref[0], 1e-3, "elastica v/L" );
    check( -tip( 0 ) / L, ref[1], 1e-3, "elastica (L-x)/L" );
    check( tip( 2 ), ref[2], 1e-3, "elastica theta" );
  }
  std::cout << "  Newton iterations: " << its << "\n";
  // Mattiasson (1981), lambda = 10: v/L 0.81061, (L-x)/L 0.55500; tip rotation 1.430286 (elliptic integral,
  // independently by scipy)
  const auto ref10 = elastica( 10 );
  check( ref10[0], 0.81061, 2e-5, "shooting vs Mattiasson v/L" );
  check( ref10[1], 0.55500, 2e-5, "shooting vs Mattiasson (L-x)/L" );
  check( ref10[2], 1.430286, 2e-6, "shooting vs elliptic integral theta" );
}

/// von Mises section points: the moment-curvature path of a beam in pure bending equals B2D2 at small rotations and
/// is independent of a superposed large rigid rotation; the plastic moment is reached
void testPlasticBending()
{
  const double                L = 100, fy = 500, hgt = 20, b = 10;
  const auto                  section = rectangle2D( hgt, b, 10 );
  const std::vector< double > vm      = { E, nu, fy, 0., 0., 1., 7.85e-9 };
  const double                Mp      = fy * b * hgt * hgt / 4;
  const double                kappaY  = 2 * fy / ( E * hgt );

  Beam            lin( "B2D2", { 0., 0., L, 0. }, section, "VONMISES", vm );
  Beam            cr( "B2D2CR", { 0., 0., L, 0. }, section, "VONMISES", vm );
  Beam            crRot( "B2D2CR", { 0., 0., L, 0. }, section, "VONMISES", vm );
  Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 6 ), URotOld = UOld;
  double          maxDiff = 0, maxDiffRot = 0;
  std::cout << "kappa/kappa_y  M_B2D2/M_p  M_B2D2CR/M_p\n";
  for ( int i = 1; i <= 40; i++ ) {
    // load to 8 kappa_y, unload to -4 kappa_y
    const double    k     = i <= 20 ? 8 * kappaY * i / 20. : 8 * kappaY - 12 * kappaY * ( i - 20 ) / 20.;
    const double    Theta = 0.5 * k * L;
    Eigen::VectorXd U( 6 );
    U << 0, 0, -Theta, 0, 0, Theta;
    lin.kernels( U, UOld, true );
    cr.kernels( U, UOld, true );
    const Eigen::VectorXd URot = rigidMotion( { 0., 0. }, { L, 0. }, U, 2.0 + 0.1 * i, { 1., 2. } );
    crRot.kernels( URot, URotOld, true );
    UOld = U, URotOld = URot;
    const double Ml = lin.state( "bending moment", 1 )( 0 ), Mc = cr.state( "bending moment", 1 )( 0 ),
                 Mr = crRot.state( "bending moment", 1 )( 0 );
    maxDiff         = std::max( maxDiff, std::abs( Ml - Mc ) / Mp );
    maxDiffRot      = std::max( maxDiffRot, std::abs( Mr - Mc ) / Mp );
    if ( i % 5 == 0 )
      std::cout << std::setprecision( 6 ) << k / kappaY << "  " << Ml / Mp << "  " << Mc / Mp << "\n";
  }
  std::cout << "max |M_CR - M_B2D2| / M_p = " << maxDiff << ", rigidly rotating CR vs CR: " << maxDiffRot << "\n";
  // pure bending with the nodes fixed: the chord does not rotate, the local rotations are exact -> identical
  throwExceptionOnFailure( maxDiff < 1e-12 && maxDiffRot < 1e-9, "plastic bending: CR = B2D2" );

  // a plastic cantilever at small rotation: tip load path; CR and B2D2 agree to the order of the rotation
  const auto            sec2 = rectangle2D( hgt, b, 5 );
  const double          Py   = fy * b * hgt * hgt / 6 / ( 2 * L ); // first yield at the root of a 2L cantilever
  std::vector< double > v, diss;
  for ( const char* type : { "B2D2", "B2D2CR" } ) {
    Chain chain( type, 8, 2 * L, sec2, "VONMISES", { E, nu, fy, 2000., 0., 1., 7.85e-9 } );
    chain.solve( { 0, 1.3 * Py, 0 }, 13, 1e-9 );
    v.push_back( chain.tip()( 1 ) );
    diss.push_back( chain.beams[0].state( "dissipation", 0 )( 0 ) );
    std::cout << std::setprecision( 8 ) << type << ": plastic cantilever 1.3 P_y: tip v = " << v.back() << ", rotation "
              << chain.tip()( 2 ) << ", root dissipation " << diss.back() << "\n";
  }
  throwExceptionOnFailure( diss[0] > 0 && std::abs( v[1] - v[0] ) < 1e-3 * v[0] &&
                             std::abs( diss[1] - diss[0] ) < 2e-2 * diss[0],
                           "plastic cantilever at small rotation: CR != B2D2" );
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactory,
    testRigidBodyRotations,
    testSmallDeformationLimit,
    testTangents,
    testRollUp,
    testElastica,
    testPlasticBending,
  } );
  return 0;
}
