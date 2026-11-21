#include "Marmot/GM7.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotMaterialGradientEnhancedCorotational.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material       = Marmot::Materials::MaterialGradientEnhancedFiniteStrain< 6 >;
  using Response       = Material::ConstitutiveResponse< 3 >;
  using Tangents       = Material::AlgorithmicModuli< 3 >;
  using Deformation    = Material::Deformation< 3 >;
  using TimeIncrement  = Material::TimeIncrement;
  using NonlocalVector = Material::NonlocalVector;
  using Tensor33d      = Marmot::Math::Tensors::FastorTensors::Tensor33d;

  // {E, nu, k1, k2, k3, k4, k5, c, m, c3, fc, transitionWidth, transitionExponent, stressUnit, density, eta}
  const std::array< double, 16 > properties =
    { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4, 25.0, 1.0, 4.0, 36.0, 4.0, 0.4, 1.0, 2.4e-9, 3e-4 };

  const std::string where = " in " + std::string( __FILE__ );

  std::unique_ptr< Material > createCorotationalGM7()
  {
    return std::unique_ptr< Material >(
      Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory< 6 >::createMaterial( "GM7/COROTATIONAL",
                                                                                          properties.data(),
                                                                                          int( properties.size() ),
                                                                                          1 ) );
  }

  Eigen::Matrix3d toEigen( const Tensor33d& t )
  {
    return Eigen::Matrix3d( Marmot::Math::Tensors::mapEigenToFastor( t ) );
  }

  Tensor33d toFastor( const Eigen::Matrix3d& m )
  {
    Tensor33d t;
    Marmot::Math::Tensors::mapEigenToFastor( t ) = m;
    return t;
  }

  /// the logarithmic strain ln U = 1/2 ln(F^T F) in Voigt notation (engineering shear)
  NonlocalVector logarithmicStrain( const Eigen::Matrix3d& F )
  {
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigenC( F.transpose() * F );
    const Eigen::Matrix3d                                  logU = eigenC.eigenvectors() *
                                 ( 0.5 * eigenC.eigenvalues().array().log() ).matrix().asDiagonal() *
                                 eigenC.eigenvectors().transpose();
    return Marmot::ContinuumMechanics::Voigt::voigtFromStrainMatrix< 3 >( logU );
  }

  Eigen::Matrix3d rotation( double angle )
  {
    return Eigen::AngleAxisd( angle, Eigen::Vector3d( 1., 2., 3. ).normalized() ).toRotationMatrix();
  }

  /// an anisotropic stretch path that runs into tensile damage
  Eigen::Matrix3d stretch( int step, int nSteps )
  {
    Eigen::Matrix3d E;
    E << 1.2e-3, 1.8e-4, 0.0, 1.8e-4, -2.4e-4, 6e-5, 0.0, 6e-5, -2.4e-4;
    return Eigen::Matrix3d::Identity() + double( step ) / nSteps * E;
  }

  struct Point {
    std::unique_ptr< Material > material = createCorotationalGM7();
    std::vector< double >       stateVars;
    Response                    response;

    Point()
    {
      stateVars.assign( material->getNumberOfRequiredStateVars(), 0.0 );
      material->initializeYourself( stateVars.data(), int( stateVars.size() ) );
      response.stateVars = stateVars.data();
    }

    // the nonlocal strain of a homogeneous field equals the local one, the logarithmic strain
    Tangents step( const Eigen::Matrix3d& F, double time )
    {
      Tangents tangents;
      material->computeStress( response,
                               tangents,
                               Deformation{ toFastor( F ), logarithmicStrain( F ) },
                               TimeIncrement{ time, 1.0 } );
      return tangents;
    }
  };

} // namespace

// The Daleckii-Krein derivative of 1/2 ln C against central differences, for distinct and for repeated eigenvalues.
void testSpectralDerivative()
{
  using Wrapper = Marmot::Materials::GradientEnhancedCorotationalWrapper< Marmot::Materials::GM7 >;

  const auto logC = []( const Eigen::Matrix3d& C ) {
    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > e( C );
    return Eigen::Matrix3d( e.eigenvectors() * ( 0.5 * e.eigenvalues().array().log() ).matrix().asDiagonal() *
                            e.eigenvectors().transpose() );
  };

  const Eigen::Matrix3d Q = rotation( 0.7 );
  for ( const Eigen::Vector3d lambda : { Eigen::Vector3d( 1.3, 0.8, 1.05 ), Eigen::Vector3d( 1.2, 1.2, 0.9 ) } ) {
    const Eigen::Matrix3d C = Q * lambda.asDiagonal() * Q.transpose();

    const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > e( C );
    const auto dLogU_dC = Wrapper::spectralDerivative( e.eigenvalues(), e.eigenvectors(), []( double a, double b ) {
      const double x = ( a - b ) / b;
      return 0.5 * ( x == 0.0 ? 1.0 / b : std::log1p( x ) / ( a - b ) );
    } );

    const double h = 1e-6;
    for ( int I = 0; I < 3; I++ )
      for ( int J = 0; J < 3; J++ ) {
        // a symmetric perturbation of C_IJ and C_JI
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
        H( I, J ) += 0.5;
        H( J, I ) += 0.5;
        const Eigen::Matrix3d numerical = ( logC( C + h * H ) - logC( C - h * H ) ) / ( 2 * h );
        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ ) {
            const double analytical = 0.5 * ( dLogU_dC( i, j, I, J ) + dLogU_dC( i, j, J, I ) );
            throwExceptionOnFailure( std::abs( analytical - numerical( i, j ) ) < 1e-8,
                                     Marmot::Utilities::MakeString()
                                       << "d(ln U)/dC(" << i << j << I << J << ") = " << analytical << " vs "
                                       << numerical( i, j ) << where );
          }
      }
  }
}

// A rigid rotation superposed on every increment of a damaging, anisotropic stretch path rotates the Kirchhoff stress
// exactly and leaves the local driving forces unchanged: the test that the wrapper exists for.
void testRigidRotationRotatesTheStressExactly()
{
  Point     reference, rotated;
  const int nSteps = 20;

  double maxStress = 0.0;
  for ( int k = 1; k <= nSteps; k++ ) {
    const Eigen::Matrix3d U = stretch( k, nSteps );
    const Eigen::Matrix3d Q = rotation( 0.07 * k );

    reference.step( U, k );
    rotated.step( Q * U, k );

    const Eigen::Matrix3d tau        = toEigen( reference.response.tau );
    const Eigen::Matrix3d tauRotated = toEigen( rotated.response.tau );
    maxStress                        = std::max( maxStress, tau( 0, 0 ) );

    throwExceptionOnFailure( ( tauRotated - Q * tau * Q.transpose() ).norm() < 1e-10 * tau.norm(),
                             Marmot::Utilities::MakeString()
                               << "step " << k << ": the stress did not rotate with the rigid rotation, error "
                               << ( tauRotated - Q * tau * Q.transpose() ).norm() << where );
    throwExceptionOnFailure( ( rotated.response.L - reference.response.L ).norm() < 1e-10 * reference.response.L.norm(),
                             "the local driving forces changed under a rigid rotation" + where );
  }

  // the path must have reached the softening branch, so that the damaged response is rotated, too
  throwExceptionOnFailure( toEigen( reference.response.tau )( 0, 0 ) < 0.95 * maxStress,
                           Marmot::Utilities::MakeString()
                             << "the stretch path did not soften: tau_xx " << toEigen( reference.response.tau )( 0, 0 )
                             << " at the end, " << maxStress << " at most" << where );
}

// Under a pure stretch the wrapped GM7 sees the logarithmic strain ln(lambda) = (lambda - 1) + O((lambda - 1)^2), so
// the Kirchhoff stress equals the stress of small-strain GM7 for the engineering strain up to second order: halving the
// strain divides the difference by four.
void testPureStretchIsSmallStrainGM7()
{
  auto differenceToSmallStrain = []( double strain ) {
    const int nSteps = 10;

    Point finiteStrain;

    std::unique_ptr< Marmot::Materials::MaterialGradientEnhancedHypoElastic< 6 > > smallStrain(
      Marmot::Registry::MaterialGradientEnhancedHypoElasticFactory< 6 >::createMaterial( "GM7",
                                                                                         properties.data(),
                                                                                         int( properties.size() ),
                                                                                         1 ) );
    std::vector< double > stateVars( smallStrain->getNumberOfRequiredStateVars(), 0.0 );
    Marmot::Materials::MaterialGradientEnhancedHypoElastic< 6 >::response  res;
    Marmot::Materials::MaterialGradientEnhancedHypoElastic< 6 >::tangents  tan;
    Marmot::Materials::MaterialGradientEnhancedHypoElastic< 6 >::increment inc;
    res.stress.setZero();
    res.stateVars = stateVars.data();

    for ( int k = 1; k <= nSteps; k++ ) {
      Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
      F( 0, 0 ) += strain * k / nSteps;
      finiteStrain.step( F, k );

      inc.dStrain      = Marmot::Math::Vector6d::Zero();
      inc.dStrain( 0 ) = strain / nSteps;
      inc.K            = Eigen::Vector< double, 6 >::Zero();
      inc.K( 0 )       = strain * k / nSteps;
      inc.dK           = inc.dStrain;
      inc.time         = k;
      inc.dT           = 1.0;
      smallStrain->computeStress( res, tan, inc );
    }

    const Eigen::Matrix3d tau   = toEigen( finiteStrain.response.tau );
    const Eigen::Matrix3d sigma = Marmot::ContinuumMechanics::Voigt::stressMatrixFromVoigt< 3 >( res.stress );
    return std::make_pair( ( tau - sigma ).norm(), sigma.norm() );
  };

  const auto [error, stress]         = differenceToSmallStrain( 1e-4 );
  const auto [errorHalf, stressHalf] = differenceToSmallStrain( 5e-5 );

  throwExceptionOnFailure( error < 1e-4 * stress, "pure stretch: not the small-strain GM7 response" + where );
  throwExceptionOnFailure( error / errorHalf > 3.5 && error / errorHalf < 4.5,
                           Marmot::Utilities::MakeString() << "pure stretch: the difference to small-strain GM7 is not "
                                                              "of second order, ratio "
                                                           << error / errorHalf << where );
}

// The analytic tangent against central differences of the stress update, at a damaged and rotated state.
void testConsistentTangent()
{
  Point     point;
  const int nSteps = 20;
  for ( int k = 1; k <= 14; k++ )
    point.step( rotation( 0.07 * k ) * stretch( k, nSteps ), k );

  const Eigen::Matrix3d F = rotation( 0.07 * 15 ) * stretch( 15, nSteps );
  const NonlocalVector  N = logarithmicStrain( F );

  const std::vector< double > stateVarsOld = point.stateVars;
  const Response              responseOld  = point.response;

  auto evaluate = [&]( const Eigen::Matrix3d& F_, const NonlocalVector& N_ ) {
    point.stateVars    = stateVarsOld;
    Response response  = responseOld;
    response.stateVars = point.stateVars.data();
    Tangents tangents;
    point.material->computeStress( response, tangents, Deformation{ toFastor( F_ ), N_ }, TimeIncrement{ 15, 1.0 } );
    return std::make_tuple( toEigen( response.tau ), NonlocalVector( response.L ), tangents );
  };

  const auto [tau, L, tangents] = evaluate( F, N );
  const double scale            = tau.norm() / 1e-4; // stiffness scale

  const double h = 1e-8;
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Eigen::Matrix3d dF        = Eigen::Matrix3d::Zero();
      dF( k, l )                = h;
      const auto [tauP, LP, tP] = evaluate( F + dF, N );
      const auto [tauM, LM, tM] = evaluate( F - dF, N );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ ) {
          const double numerical = ( tauP( i, j ) - tauM( i, j ) ) / ( 2 * h );
          throwExceptionOnFailure( std::abs( tangents.dTau_dF( i, j, k, l ) - numerical ) < 1e-5 * scale,
                                   Marmot::Utilities::MakeString()
                                     << "dTau_dF(" << i << j << k << l << ") = " << tangents.dTau_dF( i, j, k, l )
                                     << " vs " << numerical << where );
        }
      for ( int K = 0; K < 6; K++ ) {
        const double numerical = ( LP( K ) - LM( K ) ) / ( 2 * h );
        throwExceptionOnFailure( std::abs( tangents.dL_dF[K]( k, l ) - numerical ) < 1e-6,
                                 Marmot::Utilities::MakeString()
                                   << "dL_dF[" << K << "](" << k << l << ") = " << tangents.dL_dF[K]( k, l ) << " vs "
                                   << numerical << where );
      }
    }

  for ( int M = 0; M < 6; M++ ) {
    NonlocalVector dN         = NonlocalVector::Zero();
    dN( M )                   = h;
    const auto [tauP, LP, tP] = evaluate( F, N + dN );
    const auto [tauM, LM, tM] = evaluate( F, N - dN );
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ ) {
        const double numerical = ( tauP( i, j ) - tauM( i, j ) ) / ( 2 * h );
        throwExceptionOnFailure( std::abs( tangents.dTau_dN[M]( i, j ) - numerical ) < 1e-5 * scale,
                                 Marmot::Utilities::MakeString()
                                   << "dTau_dN[" << M << "](" << i << j << ") = " << tangents.dTau_dN[M]( i, j )
                                   << " vs " << numerical << where );
      }
    for ( int K = 0; K < 6; K++ ) {
      const double numerical = ( LP( K ) - LM( K ) ) / ( 2 * h );
      throwExceptionOnFailure( std::abs( tangents.dL_dN( K, M ) - numerical ) < 1e-6, "dL_dN inconsistent" + where );
    }
  }

  // the explicit update gives the same stress, without the tangent
  point.stateVars    = stateVarsOld;
  Response response  = responseOld;
  response.stateVars = point.stateVars.data();
  point.material->computeStressExplicit( response, Deformation{ toFastor( F ), N }, TimeIncrement{ 15, 1.0 } );
  throwExceptionOnFailure( ( toEigen( response.tau ) - tau ).norm() == 0.0,
                           "computeStressExplicit differs from computeStress" + where );
}

// The wrapper forwards GM7's density, nonlocal viscosity and micro-inertia, and resolves GM7's state names.
void testForwardedProperties()
{
  Point p;
  throwExceptionOnFailure( p.material->getDensity( p.stateVars.data() ) == properties[14], "density" + where );
  throwExceptionOnFailure( p.material->getNonlocalViscosity( p.stateVars.data() ) ==
                             NonlocalVector::Constant( properties[15] ),
                           "nonlocal viscosity" + where );
  throwExceptionOnFailure( p.material->getNonlocalMicroInertia( p.stateVars.data() ) == NonlocalVector::Zero(),
                           "micro-inertia" + where );
  throwExceptionOnFailure( p.material->getStateView( "strain_macro0Local", p.stateVars.data() ).stateSize == 6,
                           "the wrapped material's states must be reachable" + where );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testSpectralDerivative,
                                                       testRigidRotationRotatesTheStressExactly,
                                                       testPureStretchIsSmallStrainGM7,
                                                       testConsistentTangent,
                                                       testForwardedProperties };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
