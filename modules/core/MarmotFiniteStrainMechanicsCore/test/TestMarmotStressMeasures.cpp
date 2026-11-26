#include "Marmot/MarmotFastorTensorBasics.h" // Provides standard tensor types (e.g., Tensor33d, Tensor3333d)
#include "Marmot/MarmotStressMeasures.h"     // Contains stress measure transformations like KirchhoffStressFromPK2
#include "Marmot/MarmotTesting.h"            // Includes test utilities (e.g., assertions, test execution)
#include <Fastor/tensor/Tensor.h>            // Fastor tensor library
#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::StressMeasures;
using namespace Marmot::Math::Tensors::FastorTensors;

void testKirchhoffStressFromPK2_identityTensor()
{
  using namespace Fastor;

  // Define PK2 as the 3x3 identity tensor
  Tensor33d PK2;
  PK2.eye();
  // PK2 =
  // [ 1.0  0.0  0.0 ]
  // [ 0.0  1.0  0.0 ]
  // [ 0.0  0.0  1.0 ]

  // Define F also as the identity tensor
  Tensor33d F;
  F.eye();
  // F =
  // [ 1.0  0.0  0.0 ]
  // [ 0.0  1.0  0.0 ]
  // [ 0.0  0.0  1.0 ]

  // Compute tau = F * PK2 * F^T
  Tensor33d tau = KirchhoffStressFromPK2( PK2, F );

  // Since both F and PK2 are identity, tau should equal PK2
  Tensor33d expectedTau = PK2;
  // expectedTau =
  // [ 1.0  0.0  0.0 ]
  // [ 0.0  1.0  0.0 ]
  // [ 0.0  0.0  1.0 ]

  throwExceptionOnFailure( checkIfEqual( tau, expectedTau, 1e-12 ),
                           "KirchhoffStressFromPK2 failed while using Identity Tensor in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressFromPK2_diagonalTensor()
{
  using namespace Fastor;

  // Define a diagonal PK2 tensor
  Tensor33d PK2_dt( 0.0 );
  PK2_dt( 0, 0 ) = 2.0;
  PK2_dt( 1, 1 ) = 1.5;
  PK2_dt( 2, 2 ) = 1.0;
  // PK2_dt =
  // [ 2.0  0.0  0.0 ]
  // [ 0.0  1.5  0.0 ]
  // [ 0.0  0.0  1.0 ]

  // Define a diagonal deformation gradient F
  Tensor33d F_dt( 0.0 );
  F_dt( 0, 0 ) = 1.2;
  F_dt( 1, 1 ) = 0.9;
  F_dt( 2, 2 ) = 1.1;
  // F_dt =
  // [ 1.2  0.0  0.0 ]
  // [ 0.0  0.9  0.0 ]
  // [ 0.0  0.0  1.1 ]

  // Compute tau = F * PK2 * F^T
  Tensor33d tau_dt = KirchhoffStressFromPK2( PK2_dt, F_dt );

  // Manually compute expected result
  Tensor33d expectedTau_dt( 0.0 );
  expectedTau_dt( 0, 0 ) = 1.2 * 2.0 * 1.2; // 2.88
  expectedTau_dt( 1, 1 ) = 0.9 * 1.5 * 0.9; // 1.215
  expectedTau_dt( 2, 2 ) = 1.1 * 1.0 * 1.1; // 1.21
  // expectedTau_dt =
  // [ 2.88   0.0     0.0   ]
  // [ 0.0    1.215   0.0   ]
  // [ 0.0    0.0     1.21  ]

  throwExceptionOnFailure( checkIfEqual( tau_dt, expectedTau_dt, 1e-12 ),
                           "KirchhoffStressFromPK2 failed while using Diagonal Tensor in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressFromPK2_offDiagonalTensor()
{
  using namespace Fastor;

  // Define a symmetric PK2 tensor with off-diagonal terms
  Tensor33d PK2_off( 0.0 );
  PK2_off( 0, 0 ) = 2.0;
  PK2_off( 0, 1 ) = 0.5;
  PK2_off( 0, 2 ) = 0.1;
  PK2_off( 1, 0 ) = 0.5;
  PK2_off( 1, 1 ) = 1.5;
  PK2_off( 1, 2 ) = 0.3;
  PK2_off( 2, 0 ) = 0.1;
  PK2_off( 2, 1 ) = 0.3;
  PK2_off( 2, 2 ) = 1.0;
  // PK2_off =
  // [ 2.0  0.5  0.1 ]
  // [ 0.5  1.5  0.3 ]
  // [ 0.1  0.3  1.0 ]

  // Define a general deformation gradient F with off-diagonal entries
  Tensor33d F_off( 0.0 );
  F_off( 0, 0 ) = 1.1;
  F_off( 0, 1 ) = 0.2;
  F_off( 0, 2 ) = 0.0;
  F_off( 1, 0 ) = 0.0;
  F_off( 1, 1 ) = 0.9;
  F_off( 1, 2 ) = 0.1;
  F_off( 2, 0 ) = 0.0;
  F_off( 2, 1 ) = 0.0;
  F_off( 2, 2 ) = 1.2;
  // F_off =
  // [ 1.1  0.2  0.0 ]
  // [ 0.0  0.9  0.1 ]
  // [ 0.0  0.0  1.2 ]

  // Compute tau = F * PK2 * F^T
  Tensor33d tau_off = KirchhoffStressFromPK2( PK2_off, F_off );

  // Calculate expected tau manually
  Tensor33d expectedTau_off( 0.0 );
  expectedTau_off( 0, 0 ) = 2.70;
  expectedTau_off( 0, 1 ) = 0.782;
  expectedTau_off( 0, 2 ) = 0.204;
  expectedTau_off( 1, 0 ) = 0.782;
  expectedTau_off( 1, 1 ) = 1.279;
  expectedTau_off( 1, 2 ) = 0.444;
  expectedTau_off( 2, 0 ) = 0.204;
  expectedTau_off( 2, 1 ) = 0.444;
  expectedTau_off( 2, 2 ) = 1.44;
  // expectedTau_off =
  // [ 2.70   0.782  0.204 ]
  // [ 0.782  1.279  0.444 ]
  // [ 0.204  0.444  1.44  ]

  throwExceptionOnFailure( checkIfEqual( tau_off, expectedTau_off, 1e-12 ),
                           "KirchhoffStressFromPK2 failed while using Off-Diagonal Tensor in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressFromPK2FirstOrderDerived_diagonalTensor()
{
  using namespace Fastor;

  // Define a diagonal PK2 tensor
  Tensor33d PK2( 0.0 );
  ;
  PK2( 0, 0 ) = 2.0;
  PK2( 1, 1 ) = 1.5;
  PK2( 2, 2 ) = 1.0;
  // PK2 =
  // [ 2.0  0.0  0.0 ]
  // [ 0.0  1.5  0.0 ]
  // [ 0.0  0.0  1.0 ]

  // Define a diagonal deformation gradient F
  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.2;
  F( 1, 1 ) = 0.9;
  F( 2, 2 ) = 1.1;
  // F =
  // [ 1.2  0.0  0.0 ]
  // [ 0.0  0.9  0.0 ]
  // [ 0.0  0.0  1.1 ]

  // Call function to compute tau, and its first-order derivatives
  auto [tau, dTau_dPK2, dTau_dF] = FirstOrder::KirchhoffStressFromPK2( PK2, F );

  // Test 1: Check tau result
  Tensor33d expectedTau( 0.0 );
  for ( int i = 0; i < 3; ++i )
    expectedTau( i, i ) = F( i, i ) * PK2( i, i ) * F( i, i );
  // expectedTau =
  // [ 2.88   0.0    0.0  ]
  // [ 0.0    1.215  0.0  ]
  // [ 0.0    0.0    1.21 ]

  throwExceptionOnFailure( checkIfEqual( tau, expectedTau, 1e-12 ),
                           "FirstOrder tau with Diagonal Tensor mismatch in " + std::string( __PRETTY_FUNCTION__ ) );

  // Test 2: Check dTau/dPK2
  Tensor3333d expected_dTau_dPK2( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q )
          expected_dTau_dPK2( i, j, p, q ) = F( i, p ) * F( j, q );

  throwExceptionOnFailure( checkIfEqual( dTau_dPK2, expected_dTau_dPK2, 1e-12 ),
                           "FirstOrder dTau_dPK2 with Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Test 3: Check dTau/dF
  Tensor3333d expected_dTau_dF( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q ) {
          double term1                   = ( i == p ) ? PK2( q, q ) * F( j, q ) : 0.0;
          double term2                   = ( j == p ) ? F( i, q ) * PK2( q, q ) : 0.0;
          expected_dTau_dF( i, j, p, q ) = term1 + term2;
        }

  throwExceptionOnFailure( checkIfEqual( dTau_dF, expected_dTau_dF, 1e-12 ),
                           "FirstOrder dTau_dF with Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressFromPK2FirstOrderDerived_offDiagonalTensor()
{
  using namespace Fastor;

  // Define symmetric PK2 tensor
  Tensor33d PK2( 0.0 );
  PK2( 0, 0 ) = 2.0;
  PK2( 0, 1 ) = 1.0;
  PK2( 0, 2 ) = 2.0;
  PK2( 1, 0 ) = 1.0;
  PK2( 1, 1 ) = 3.0;
  PK2( 1, 2 ) = 2.0;
  PK2( 2, 0 ) = 2.0;
  PK2( 2, 1 ) = 2.0;
  PK2( 2, 2 ) = 1.0;
  // PK2 =
  // [ 2.0  1.0  2.0 ]
  // [ 1.0  3.0  2.0 ]
  // [ 2.0  2.0  1.0 ]

  // Define general F
  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.0;
  F( 0, 1 ) = 2.0;
  F( 0, 2 ) = 0.0;
  F( 1, 0 ) = 0.0;
  F( 1, 1 ) = 3.0;
  F( 1, 2 ) = 1.0;
  F( 2, 0 ) = 0.0;
  F( 2, 1 ) = 0.0;
  F( 2, 2 ) = 2.0;
  // F =
  // [ 1.0  2.0  0.0 ]
  // [ 0.0  3.0  1.0 ]
  // [ 0.0  0.0  2.0 ]

  // Call analytical function
  auto [tau, dTau_dPK2, dTau_dF] = FirstOrder::KirchhoffStressFromPK2( PK2, F );

  // Test 1: Validate tau
  auto computeTau = [&]( const Tensor33d& A, const Tensor33d& B ) {
    Tensor33d result( 0.0 );
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l )
            result( i, j ) += A( i, k ) * B( k, l ) * A( j, l );
    return result;
  };

  Tensor33d tau_ref = computeTau( F, PK2 );

  throwExceptionOnFailure( checkIfEqual( tau, tau_ref, 1e-12 ),
                           "FirstOrder tau with Off-Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Test 2: Validate dTau/dPK2 analytically
  Tensor3333d ref_dTau_dPK2( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q )
          ref_dTau_dPK2( i, j, p, q ) = F( i, p ) * F( j, q );

  throwExceptionOnFailure( checkIfEqual( dTau_dPK2, ref_dTau_dPK2, 1e-12 ),
                           "FirstOrder dTau_dPK2 with Off-Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Test 3: Check dTau/dF analytically
  Tensor3333d expected_dTau_dF( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q ) {
          double sum1 = 0.0;
          double sum2 = 0.0;
          if ( i == p ) {
            for ( int k = 0; k < 3; ++k )
              sum1 += PK2( q, k ) * F( j, k );
          }
          if ( j == p ) {
            for ( int k = 0; k < 3; ++k )
              sum2 += F( i, k ) * PK2( k, q );
          }
          expected_dTau_dF( i, j, p, q ) = sum1 + sum2;
        }

  throwExceptionOnFailure( checkIfEqual( dTau_dF, expected_dTau_dF, 1e-12 ),
                           "FirstOrder dTau_dF with Off-Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

namespace {
  /* Deterministic, non-symmetric test data. */
  Tensor33d testDeformationGradient()
  {
    return Tensor33d{ { 1.20, 0.10, -0.03 }, { 0.05, 0.90, 0.07 }, { -0.02, 0.04, 1.05 } };
  }

  Tensor33d testPK2()
  {
    return Tensor33d{ { 3.0, 0.4, -0.7 }, { 0.9, -2.0, 0.3 }, { 0.2, -0.5, 1.5 } };
  }

  Tensor3333d testdPK2_dC()
  {
    Tensor3333d D;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l )
            D( i, j, k, l ) = std::sin( 1. + 27. * i + 9. * j + 3. * k + l ) + 0.1 * ( i - 2. * j + 3. * k - l );
    return D;
  }
} // namespace

void testKirchhoffStressAndTangentFromPK2_matchesChainOfDerivatives()
{
  using namespace Fastor;
  using namespace Marmot::Math::Tensors::FastorIndices;

  const Tensor33d   F = testDeformationGradient();
  const Tensor3333d D = testdPK2_dC();
  // the partial derivative dTau/dF of FirstOrder::KirchhoffStressFromPK2 holds for a symmetric PK2
  const Tensor33d S = 0.5 * ( testPK2() + transpose( testPK2() ) );

  const auto [tauRef, dTau_dPK2, dTau_dFPartial] = FirstOrder::KirchhoffStressFromPK2( S, F );
  const Tensor3333d dC_dF = einsum< IK, kJ, to_IJkK >( Spatial3D::I, F ) + einsum< kI, JK, to_IJkK >( F, Spatial3D::I );
  const Tensor3333d dTau_dFRef = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, D ), dC_dF ) + dTau_dFPartial;

  const auto [tau, dTau_dF] = FirstOrder::KirchhoffStressAndTangentFromPK2( S, D, F );

  throwExceptionOnFailure( checkIfEqual( tau, tauRef, 1e-13 ),
                           "tau mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( dTau_dF, dTau_dFRef, 1e-12 ),
                           "dTau_dF mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressAndTangentFromPK2_matchesFiniteDifferences()
{
  // S( C ) = S0 + D : ( C - C0 ) with non-symmetric S0 and D: the tangent is the total derivative of
  // tau( F ) = F S( F^T F ) F^T, checked by central differences
  const Tensor33d   F0 = testDeformationGradient();
  const Tensor33d   S0 = testPK2();
  const Tensor3333d D  = testdPK2_dC();
  const Tensor33d   C0 = Fastor::matmul( Fastor::transpose( F0 ), F0 );

  const auto tauOf = [&]( const Tensor33d& F ) {
    const Tensor33d dC = Fastor::matmul( Fastor::transpose( F ), F ) - C0;
    Tensor33d       S  = S0;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l )
            S( i, j ) += D( i, j, k, l ) * dC( k, l );
    return KirchhoffStressFromPK2( S, F );
  };

  const auto [tau, dTau_dF] = FirstOrder::KirchhoffStressAndTangentFromPK2( S0, D, F0 );

  throwExceptionOnFailure( checkIfEqual( tau, tauOf( F0 ), 1e-13 ),
                           "tau mismatch in " + std::string( __PRETTY_FUNCTION__ ) );

  const double h = 1e-6;
  Tensor3333d  dTau_dFNumeric;
  for ( int k = 0; k < 3; ++k )
    for ( int M = 0; M < 3; ++M ) {
      Tensor33d Fp = F0, Fm = F0;
      Fp( k, M ) += h;
      Fm( k, M ) -= h;
      const Tensor33d dTau = ( tauOf( Fp ) - tauOf( Fm ) ) / ( 2. * h );
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          dTau_dFNumeric( i, j, k, M ) = dTau( i, j );
    }

  throwExceptionOnFailure( checkIfEqual( dTau_dF, dTau_dFNumeric, 1e-7 ),
                           "dTau_dF mismatch in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  // Register all test cases
  auto tests = std::vector< std::function< void() > >{ testKirchhoffStressFromPK2_identityTensor,
                                                       testKirchhoffStressFromPK2_diagonalTensor,
                                                       testKirchhoffStressFromPK2_offDiagonalTensor,
                                                       testKirchhoffStressFromPK2FirstOrderDerived_diagonalTensor,
                                                       testKirchhoffStressFromPK2FirstOrderDerived_offDiagonalTensor,
                                                       testKirchhoffStressAndTangentFromPK2_matchesChainOfDerivatives,
                                                       testKirchhoffStressAndTangentFromPK2_matchesFiniteDifferences };

  // Execute and report any exceptions from failed tests
  executeTestsAndCollectExceptions( tests );
  return 0;
}
