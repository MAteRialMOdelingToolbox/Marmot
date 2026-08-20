#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/Ogden.h"
#include <cmath>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

namespace {

  // builds materialProperties = [nTerms, mu_1..mu_N, alpha_1..alpha_N]
  std::vector< double > makeProperties( const std::vector< double >& mu, const std::vector< double >& alpha )
  {
    std::vector< double > props;
    props.push_back( static_cast< double >( mu.size() ) );
    props.insert( props.end(), mu.begin(), mu.end() );
    props.insert( props.end(), alpha.begin(), alpha.end() );
    return props;
  }

  Ogden::ConstitutiveResponse< 3 > computeResponse( const std::vector< double >&   mu,
                                                    const std::vector< double >&   alpha,
                                                    const Tensor33d&               F,
                                                    Ogden::AlgorithmicModuli< 3 >* tangentOut = nullptr )
  {
    std::vector< double > props = makeProperties( mu, alpha );
    const Ogden           mat( props.data(), static_cast< int >( props.size() ), 1 );

    Ogden::Deformation< 3 >          def     = { F };
    Ogden::TimeIncrement             timeInc = { 0, 0.1 };
    Ogden::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
    Ogden::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };

    mat.computeStress( response, tangent, def, timeInc );

    if ( tangentOut )
      *tangentOut = tangent;

    return response;
  }

} // namespace

// Test O-1: at the undeformed configuration, the isochoric Ogden stress must vanish for any
// choice of parameters (lambda_bar_i = 1 for all i).
void testUndeformedResponse()
{
  const std::vector< double > mu    = { 1500., 300., -50. };
  const std::vector< double > alpha = { 2., 4., -2. };

  const Tensor33d inputF   = Spatial3D::I;
  const auto      response = computeResponse( mu, alpha, inputF );

  throwExceptionOnFailure( checkIfEqual( response.tau, Tensor33d( 0.0 ), 1e-10 ),
                           "O-1: Undeformed configuration must yield zero stress for Ogden material" );
}

// Test O-2: classical closed-form solution for incompressible uniaxial tension (e.g. Treloar /
// Holzapfel, "Nonlinear Solid Mechanics", Eq. 6.98/6.100). For F = diag(lambda, lambda^-1/2,
// lambda^-1/2) (J=1 exactly), the principal Kirchhoff stresses are
//   tau_1 = 2/3 * (A - B),   tau_2 = tau_3 = -1/3 * (A - B),
// with A(lambda) = sum_p mu_p lambda^alpha_p, B(lambda) = sum_p mu_p lambda^(-alpha_p/2).
void testIncompressibleUniaxialTension()
{
  const std::vector< double > mu    = { 1500., 300. };
  const std::vector< double > alpha = { 2., -2. };

  for ( const double lambda : { 0.7, 0.9, 1.0, 1.2, 1.5 } ) {

    Tensor33d inputF( 0.0 );
    inputF( 0, 0 ) = lambda;
    inputF( 1, 1 ) = 1. / std::sqrt( lambda );
    inputF( 2, 2 ) = 1. / std::sqrt( lambda );

    double A = 0., B = 0.;
    for ( size_t p = 0; p < mu.size(); ++p ) {
      A += mu[p] * std::pow( lambda, alpha[p] );
      B += mu[p] * std::pow( lambda, -alpha[p] / 2. );
    }

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 2. / 3. * ( A - B );
    stressTarget( 1, 1 ) = -1. / 3. * ( A - B );
    stressTarget( 2, 2 ) = -1. / 3. * ( A - B );

    const auto response = computeResponse( mu, alpha, inputF );

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-8 ),
                             "O-2: Incompressible uniaxial tension closed-form check failed for lambda=" +
                               std::to_string( lambda ) );

    throwExceptionOnFailure( std::abs( Fastor::trace( response.tau ) ) < 1e-8,
                             "O-2: Ogden Kirchhoff stress must be trace-free (purely isochoric) for lambda=" +
                               std::to_string( lambda ) );
  }
}

// Test O-3: algorithmic tangent dTau/dF checked against a central-difference approximation.
void testAlgorithmicTangent()
{
  const std::vector< double > mu    = { 1500., 300. };
  const std::vector< double > alpha = { 2., -2. };

  Tensor33d inputF( 0.0 );
  inputF( 0, 0 ) = 1.01;
  inputF( 0, 1 ) = 0.06;
  inputF( 0, 2 ) = -0.03;
  inputF( 1, 0 ) = 0.06;
  inputF( 1, 1 ) = 1.02;
  inputF( 1, 2 ) = 0.04;
  inputF( 2, 0 ) = -0.03;
  inputF( 2, 1 ) = 0.04;
  inputF( 2, 2 ) = 0.95;

  Ogden::AlgorithmicModuli< 3 > tangent  = { Tensor3333d( 0.0 ) };
  const auto                    response = computeResponse( mu, alpha, inputF, &tangent );

  const double h = 1e-6;
  for ( int k = 0; k < 3; ++k ) {
    for ( int K = 0; K < 3; ++K ) {
      Tensor33d Fp = inputF;
      Tensor33d Fm = inputF;
      Fp( k, K ) += h;
      Fm( k, K ) -= h;

      const auto responseP = computeResponse( mu, alpha, Fp );
      const auto responseM = computeResponse( mu, alpha, Fm );

      for ( int i = 0; i < 3; ++i ) {
        for ( int I = 0; I < 3; ++I ) {
          const double dTau_dF_numeric = ( responseP.tau( i, I ) - responseM.tau( i, I ) ) / ( 2 * h );

          throwExceptionOnFailure( std::abs( dTau_dF_numeric - tangent.dTau_dF( i, I, k, K ) ) < 1e-4,
                                   "O-3: Algorithmic tangent mismatch at (i,I,k,K)=(" + std::to_string( i ) + "," +
                                     std::to_string( I ) + "," + std::to_string( k ) + "," + std::to_string( K ) +
                                     ")" );
        }
      }
    }
  }
  (void)response;
}

// Test O-4: objectivity (Q*F) and isotropy (F*Q) for a purely isochoric arbitrary deformation.
void testRotation()
{
  const std::vector< double > mu    = { 1500., 300. };
  const std::vector< double > alpha = { 2., -2. };

  Tensor33d inputF( 0.0 );
  inputF( 0, 0 ) = 1.01;
  inputF( 0, 1 ) = 0.06;
  inputF( 0, 2 ) = -0.03;
  inputF( 1, 0 ) = 0.06;
  inputF( 1, 1 ) = 1.02;
  inputF( 1, 2 ) = 0.04;
  inputF( 2, 0 ) = -0.03;
  inputF( 2, 1 ) = 0.04;
  inputF( 2, 2 ) = 0.95;

  const auto responseUnrotated = computeResponse( mu, alpha, inputF );

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1;

    // Objectivity: tau(Q*F) = Q * tau(F) * Q^T
    {
      const Tensor33d F_rotated     = einsum< ik, kj, to_ij >( Q, inputF );
      const auto      response      = computeResponse( mu, alpha, F_rotated );
      const Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, responseUnrotated.tau, Q );

      throwExceptionOnFailure( checkIfEqual( response.tau, stressRotated, 1e-8 ),
                               "O-4a: Objectivity test failed (phi_deg=" + std::to_string( phi_deg ) + ")" );
    }

    // Isotropy: tau(F*Q) = tau(F)
    {
      const Tensor33d F_rotated = einsum< ik, kj, to_ij >( inputF, Q );
      const auto      response  = computeResponse( mu, alpha, F_rotated );

      throwExceptionOnFailure( checkIfEqual( response.tau, responseUnrotated.tau, 1e-8 ),
                               "O-4b: Isotropy test failed (phi_deg=" + std::to_string( phi_deg ) + ")" );
    }
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testUndeformedResponse,
                                                       testIncompressibleUniaxialTension,
                                                       testAlgorithmicTangent,
                                                       testRotation };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
