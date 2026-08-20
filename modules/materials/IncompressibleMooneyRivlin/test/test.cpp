#include "Marmot/IncompressibleMooneyRivlin.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

namespace {

  const double C1 = 1000.;
  const double C2 = 300.;

  IncompressibleMooneyRivlin::ConstitutiveResponse< 3 > computeResponse(
    const Tensor33d&                                    F,
    IncompressibleMooneyRivlin::AlgorithmicModuli< 3 >* tangentOut = nullptr )
  {
    const std::array< double, 2 >    props = { C1, C2 };
    const IncompressibleMooneyRivlin mat( props.data(), props.size(), 1 );

    IncompressibleMooneyRivlin::Deformation< 3 >          def     = { F };
    IncompressibleMooneyRivlin::TimeIncrement             timeInc = { 0, 0.1 };
    IncompressibleMooneyRivlin::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, nullptr );
    IncompressibleMooneyRivlin::AlgorithmicModuli< 3 >    tangent = { Tensor3333d( 0.0 ) };

    mat.computeStress( response, tangent, def, timeInc );

    if ( tangentOut )
      *tangentOut = tangent;

    return response;
  }

} // namespace

// Test M-1: at the undeformed configuration, the isochoric stress must vanish.
void testUndeformedResponse()
{
  const auto response = computeResponse( Spatial3D::I );
  throwExceptionOnFailure( checkIfEqual( response.tau, Tensor33d( 0.0 ), 1e-10 ),
                           "M-1: Undeformed configuration must yield zero stress" );
}

// Test M-2: classical closed-form incompressible uniaxial tension for Mooney-Rivlin
// (e.g. Ogden, "Non-Linear Elastic Deformations", Ch. 4). For
// F = diag(lambda, lambda^-1/2, lambda^-1/2) (J=1),
//   A(lambda) = 2*C1*lambda^2 - 2*C2*lambda^-2,
//   B(lambda) = 2*C1*lambda^-1 - 2*C2*lambda,
//   tau_1 = 2/3*(A-B), tau_2 = tau_3 = -1/3*(A-B).
void testIncompressibleUniaxialTension()
{
  for ( const double lambda : { 0.7, 0.9, 1.0, 1.2, 1.5 } ) {

    Tensor33d inputF( 0.0 );
    inputF( 0, 0 ) = lambda;
    inputF( 1, 1 ) = 1. / std::sqrt( lambda );
    inputF( 2, 2 ) = 1. / std::sqrt( lambda );

    const double A = 2. * C1 * lambda * lambda - 2. * C2 / ( lambda * lambda );
    const double B = 2. * C1 / lambda - 2. * C2 * lambda;

    Tensor33d stressTarget( 0.0 );
    stressTarget( 0, 0 ) = 2. / 3. * ( A - B );
    stressTarget( 1, 1 ) = -1. / 3. * ( A - B );
    stressTarget( 2, 2 ) = -1. / 3. * ( A - B );

    const auto response = computeResponse( inputF );

    throwExceptionOnFailure( checkIfEqual( response.tau, stressTarget, 1e-8 ),
                             "M-2: Incompressible uniaxial tension closed-form check failed for lambda=" +
                               std::to_string( lambda ) );

    throwExceptionOnFailure( std::abs( Fastor::trace( response.tau ) ) < 1e-8,
                             "M-2: Kirchhoff stress must be trace-free (purely isochoric) for lambda=" +
                               std::to_string( lambda ) );
  }
}

// Test M-3: algorithmic tangent dTau/dF checked against a central-difference approximation.
void testAlgorithmicTangent()
{
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

  IncompressibleMooneyRivlin::AlgorithmicModuli< 3 > tangent = { Tensor3333d( 0.0 ) };
  computeResponse( inputF, &tangent );

  const double h = 1e-6;
  for ( int k = 0; k < 3; ++k ) {
    for ( int K = 0; K < 3; ++K ) {
      Tensor33d Fp = inputF;
      Tensor33d Fm = inputF;
      Fp( k, K ) += h;
      Fm( k, K ) -= h;

      const auto responseP = computeResponse( Fp );
      const auto responseM = computeResponse( Fm );

      for ( int i = 0; i < 3; ++i ) {
        for ( int I = 0; I < 3; ++I ) {
          const double dTau_dF_numeric = ( responseP.tau( i, I ) - responseM.tau( i, I ) ) / ( 2 * h );

          throwExceptionOnFailure( std::abs( dTau_dF_numeric - tangent.dTau_dF( i, I, k, K ) ) < 1e-4,
                                   "M-3: Algorithmic tangent mismatch at (i,I,k,K)=(" + std::to_string( i ) + "," +
                                     std::to_string( I ) + "," + std::to_string( k ) + "," + std::to_string( K ) +
                                     ")" );
        }
      }
    }
  }
}

// Test M-4: objectivity (Q*F) and isotropy (F*Q).
void testRotation()
{
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

  const auto responseUnrotated = computeResponse( inputF );

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1;

    {
      const Tensor33d F_rotated     = einsum< ik, kj, to_ij >( Q, inputF );
      const auto      response      = computeResponse( F_rotated );
      const Tensor33d stressRotated = einsum< iI, IJ, jJ, to_ij >( Q, responseUnrotated.tau, Q );

      throwExceptionOnFailure( checkIfEqual( response.tau, stressRotated, 1e-8 ),
                               "M-4a: Objectivity test failed (phi_deg=" + std::to_string( phi_deg ) + ")" );
    }

    {
      const Tensor33d F_rotated = einsum< ik, kj, to_ij >( inputF, Q );
      const auto      response  = computeResponse( F_rotated );

      throwExceptionOnFailure( checkIfEqual( response.tau, responseUnrotated.tau, 1e-8 ),
                               "M-4b: Isotropy test failed (phi_deg=" + std::to_string( phi_deg ) + ")" );
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
