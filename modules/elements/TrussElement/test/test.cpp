#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/TrussElement.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Math;
using Marmot::Discretizations::ElementProperties;

namespace {

  /// A truss element created through the factory, together with the storage it needs
  struct Truss {
    std::unique_ptr< Marmot::Discretizations::Elements::Element > element;
    std::vector< double >                                         coordinates;
    std::vector< double >                                         properties;
    std::vector< double >                                         materialProperties;
    std::vector< double >                                         stateVars;
    int                                                           nDof;

    Truss( const std::string&    type,
           std::vector< double > coordinates_,
           double                area,
           const std::string&    material,
           std::vector< double > materialProperties_ )
      : element( Marmot::Registry::ElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( { area } ),
        materialProperties( std::move( materialProperties_ ) )
    {
      nDof = element->getNDofPerElement();
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty( ElementProperties( properties.data(), 1 ) );
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

    /// residual and tangent for the total displacement U from the committed state; commits if requested
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& U,
                                                           const Eigen::VectorXd& UOld,
                                                           bool                   commit )
    {
      const auto      backup = stateVars;
      Eigen::VectorXd P      = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K      = Eigen::MatrixXd::Zero( nDof, nDof );
      Eigen::VectorXd dU     = U - UOld;
      element->computeKernels( U.data(), dU.data(), P.data(), K.data(), 0.0, 1.0 );
      // Marmot fills K column major w.r.t. its own (row = equation) convention; it is symmetric here anyway
      if ( !commit )
        stateVars = backup;
      return { P, K };
    }

    /// central difference tangent from the committed state
    Eigen::MatrixXd numericalTangent( const Eigen::VectorXd& U, const Eigen::VectorXd& UOld, double h = 1e-7 )
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

    double state( const std::string& name, int qp = 0 ) { return *element->getStateView( name, qp ).stateLocation; }
  };

  void checkTangent( Truss&                 truss,
                     const Eigen::VectorXd& U,
                     const Eigen::VectorXd& UOld,
                     double                 relTol,
                     const std::string&     what )
  {
    const auto [P, K]  = truss.kernels( U, UOld, false );
    const auto   numK  = truss.numericalTangent( U, UOld );
    const double scale = std::max( 1.0, numK.cwiseAbs().maxCoeff() );
    throwExceptionOnFailure( checkIfEqual< double >( K, numK, relTol * scale ),
                             Marmot::Utilities::MakeString()
                               << what << ": tangent does not match the numerical one\nK=\n"
                               << K << "\nnumK=\n"
                               << numK );
  }

  const std::vector< double > linearElastic = { 210000., 0.3, 7.85e-9 };

} // namespace

/// small strain, inclined linear truss: the classic stiffness matrix EA/L [c c^T, -c c^T; ...]
void testSmallStrainInclinedLinearTruss()
{
  const double E = 210000., A = 50.;
  Truss        truss( "T2D2", { 0., 0., 300., 400. }, A, "LINEARELASTIC", linearElastic );
  const double L = 500.;

  Eigen::VectorXd U   = Eigen::VectorXd::Zero( 4 );
  U.segment< 2 >( 2 ) = Eigen::Vector2d( 0.6, 0.8 ) * 0.25; // axial elongation 0.25
  auto [P, K]         = truss.kernels( U, Eigen::VectorXd::Zero( 4 ), true );

  const Eigen::Vector2d c( 0.6, 0.8 );
  const Eigen::Matrix2d cc = c * c.transpose();
  Eigen::Matrix4d       Kexpected;
  Kexpected.block< 2, 2 >( 0, 0 ) = cc;
  Kexpected.block< 2, 2 >( 0, 2 ) = -cc;
  Kexpected.block< 2, 2 >( 2, 0 ) = -cc;
  Kexpected.block< 2, 2 >( 2, 2 ) = cc;
  Kexpected *= E * A / L;

  throwExceptionOnFailure( checkIfEqual< double >( K, Kexpected, 1e-8 * E * A / L ), "inclined truss stiffness" );
  throwExceptionOnFailure( checkIfEqual( truss.state( "normal force" ), E * A * 0.25 / L, 1e-8 ), "normal force" );
  throwExceptionOnFailure( checkIfEqual( truss.state( "strain" ), 0.25 / L, 1e-14 ), "axial strain" );
  Eigen::Vector4d Pexpected;
  Pexpected.head< 2 >() = -c * E * A * 0.25 / L;
  Pexpected.tail< 2 >() = c * E * A * 0.25 / L;
  throwExceptionOnFailure( checkIfEqual< double >( P, Pexpected, 1e-8 ), "internal forces" );
}

/// small strain, quadratic 3D truss in uniform stretch: constant strain and stress, consistent tangent, mass
void testSmallStrainQuadraticTruss3D()
{
  const double E = 210000., A = 10.;
  // node order: end, end, mid
  Truss        truss( "T3D3", { 0., 0., 0., 2., 2., 1., 1., 1., 0.5 }, A, "LINEARELASTIC", linearElastic );
  const double L = 3.;

  const Eigen::Vector3d t = Eigen::Vector3d( 2., 2., 1. ) / L;
  Eigen::VectorXd       U( 9 );
  U.segment< 3 >( 0 ) = 0. * t;
  U.segment< 3 >( 3 ) = 1e-3 * L * t;
  U.segment< 3 >( 6 ) = 0.5e-3 * L * t;
  truss.kernels( U, Eigen::VectorXd::Zero( 9 ), true );
  for ( int qp = 0; qp < 2; qp++ ) {
    throwExceptionOnFailure( checkIfEqual( truss.state( "strain", qp ), 1e-3, 1e-12 ), "uniform strain" );
    throwExceptionOnFailure( checkIfEqual( truss.state( "stress", qp ), E * 1e-3, 1e-8 ), "uniform stress" );
  }
  checkTangent( truss, U, Eigen::VectorXd::Zero( 9 ), 1e-6, "T3D3 small strain" );

  Eigen::VectorXd m = Eigen::VectorXd::Zero( 9 );
  truss.element->computeLumpedInertia( m.data() );
  throwExceptionOnFailure( checkIfEqual( m.sum() / 3, 7.85e-9 * A * L, 1e-20 ), "lumped mass" );
  throwExceptionOnFailure( checkIfEqual( m( 6 ), 2. / 3 * 7.85e-9 * A * L, 1e-20 ), "midside node mass" );
}

/// small strain with plasticity: the uniaxial reduction of a 3D von Mises model yields at the yield stress
void testSmallStrainVonMises()
{
  //                              E        nu   fy    H   dfy  delta rho
  const std::vector< double > vm = { 210000., 0.3, 500., 0., 0., 1., 7.85e-9 };
  Truss                       truss( "T2D2", { 0., 0., 100., 0. }, 1.0, "VONMISES", vm );
  Eigen::VectorXd             UOld = Eigen::VectorXd::Zero( 4 );
  for ( int i = 1; i <= 10; i++ ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( 4 );
    U( 2 )            = i * 0.05; // up to 5e-3 strain, 2x yield strain
    truss.kernels( U, UOld, true );
    UOld = U;
  }
  throwExceptionOnFailure( checkIfEqual( truss.state( "stress" ), 500., 1e-6 ), "von Mises yield in uniaxial stress" );
}

/// finite strain: large rigid body rotations are stress free, and a stretch gives the uniaxial Neo-Hooke response
void testFiniteStrainRigidRotationAndStretch()
{
  const double K = 3000., G = 1000.;
  Truss        truss( "T3D2CR", { 0., 0., 0., 1., 0., 0. }, 2.0, "COMPRESSIBLENEOHOOKE", { K, G, 1e-9 } );

  // rotate by 90 degrees about z: node 2 goes from (1,0,0) to (0,1,0)
  Eigen::VectorXd U   = Eigen::VectorXd::Zero( 6 );
  U.segment< 3 >( 3 ) = Eigen::Vector3d( -1., 1., 0. );
  auto [P, Kt]        = truss.kernels( U, Eigen::VectorXd::Zero( 6 ), false );
  throwExceptionOnFailure( P.norm() < 1e-10,
                           Marmot::Utilities::MakeString() << "rigid rotation must be stress free, |P|=" << P.norm() );

  // ... rotated and stretched by 1.5: axial force along the current axis; the lateral stresses vanish
  U.segment< 3 >( 3 ) = Eigen::Vector3d( -1., 1.5, 0. );
  auto [P2, K2]       = truss.kernels( U, Eigen::VectorXd::Zero( 6 ), true );
  const double N      = truss.state( "normal force" );
  throwExceptionOnFailure( N > 0, "tension" );
  throwExceptionOnFailure( checkIfEqual< double >( P2.segment( 3, 3 ), Eigen::Vector3d( 0., N, 0. ), 1e-9 * N ),
                           "end force along the current axis" );
  throwExceptionOnFailure( checkIfEqual< double >( P2.segment( 0, 3 ), -P2.segment( 3, 3 ), 1e-12 ), "equilibrium" );

  // Cauchy stress x current area = normal force
  const double lat = truss.state( "lateral stretches" );
  throwExceptionOnFailure( checkIfEqual( truss.state( "stress" ) * 2.0 * lat * lat, N, 1e-9 * N ),
                           "normal force = Cauchy stress x current area" );
  throwExceptionOnFailure( checkIfEqual( truss.state( "strain" ), std::log( 1.5 ), 1e-12 ), "logarithmic strain" );
}

/// finite strain consistent tangent (material and geometric part), linear and quadratic, 2D and 3D
void testFiniteStrainTangents()
{
  {
    Truss           truss( "T2D2CR", { 0., 0., 2., 1. }, 3.0, "COMPRESSIBLENEOHOOKE", { 3000., 1000., 1e-9 } );
    Eigen::VectorXd U( 4 );
    U << 0.1, -0.2, -0.5, 0.9;
    checkTangent( truss, U, Eigen::VectorXd::Zero( 4 ), 1e-6, "T2D2CR Neo-Hooke" );
  }
  {
    Truss           truss( "T3D3CR",
                           { 0., 0., 0., 2., 0., 0., 1., 0.2, 0. },
                 3.0,
                 "COMPRESSIBLENEOHOOKE",
                           { 3000., 1000., 1e-9 } );
    Eigen::VectorXd U( 9 );
    U << 0.1, -0.2, 0.05, 0.3, 0.6, -0.4, 0.2, 0.1, 0.3;
    checkTangent( truss, U, Eigen::VectorXd::Zero( 9 ), 1e-6, "T3D3CR Neo-Hooke" );
  }
}

/// finite strain plasticity: the elastoplastic tangent is consistent, and the stress is the yield stress
void testFiniteStrainJ2Plasticity()
{
  const double fy = 300.;
  //                                         K        G       fy  fyInf eta  H   impl rho
  const std::vector< double > j2 = { 175000., 80769., fy, fy, 0.0, 0., 1, 7.8e-9 };
  Truss                       truss( "T2D3CR", { 0., 0., 10., 0., 5., 0. }, 1.0, "FINITESTRAINJ2PLASTICITY", j2 );

  Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 6 );
  Eigen::VectorXd U    = Eigen::VectorXd::Zero( 6 );
  for ( int i = 1; i <= 5; i++ ) {
    U( 2 ) = i * 0.02; // up to 1 % strain
    U( 4 ) = i * 0.01;
    truss.kernels( U, UOld, true );
    UOld = U;
  }
  throwExceptionOnFailure( checkIfEqual( truss.state( "kirchhoff stress" ), fy, 1e-6 * fy ),
                           Marmot::Utilities::MakeString()
                             << "Kirchhoff stress at yield " << truss.state( "kirchhoff stress" ) );

  Eigen::VectorXd UNew = U;
  UNew( 2 ) += 0.01;
  UNew( 4 ) += 0.005;
  UNew( 3 ) += 0.3; // with rotation
  checkTangent( truss, UNew, UOld, 1e-5, "T2D3CR J2 plasticity" );
}

void testFactoryNames()
{
  for ( const auto& name : { "T2D2", "T2D3", "T3D2", "T3D3", "T2D2CR", "T2D3CR", "T3D2CR", "T3D3CR" } ) {
    auto el = std::unique_ptr< Marmot::Discretizations::Elements::Element >(
      Marmot::Registry::ElementFactory::createElement( name, 1 ) );
    throwExceptionOnFailure( el != nullptr, Marmot::Utilities::MakeString() << name << " must be registered" );
  }
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactoryNames,
    testSmallStrainInclinedLinearTruss,
    testSmallStrainQuadraticTruss3D,
    testSmallStrainVonMises,
    testFiniteStrainRigidRotationAndStretch,
    testFiniteStrainTangents,
    testFiniteStrainJ2Plasticity,
  } );
  return 0;
}
