#include "Marmot/MarmotFastorTensorBasics.h" // Provides standard tensor types (e.g., Tensor33d, Tensor3333d)
#include "Marmot/MarmotStressMeasures.h"     // Contains stress measure transformations like KirchhoffStressFromPK2
#include "Marmot/MarmotTesting.h"            // Includes test utilities (e.g., assertions, test execution)
#include <Fastor/tensor/Tensor.h>            // Fastor tensor library

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::StressMeasures;
using namespace Marmot::FastorStandardTensors;

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
  auto [tau, dTau_dPK2, dTau_dF] = FirstOrderDerived::KirchhoffStressFromPK2( PK2, F );

  // Test 1: Check tau result
  Tensor33d expectedTau( 0.0 );
  for ( int i = 0; i < 3; ++i )
    expectedTau( i, i ) = F( i, i ) * PK2( i, i ) * F( i, i );
  // expectedTau =
  // [ 2.88   0.0    0.0  ]
  // [ 0.0    1.215  0.0  ]
  // [ 0.0    0.0    1.21 ]

  throwExceptionOnFailure( checkIfEqual( tau, expectedTau, 1e-12 ),
                           "FirstOrderDerived tau with Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Test 2: Check dTau/dPK2
  Tensor3333d expected_dTau_dPK2( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q )
          expected_dTau_dPK2( i, j, p, q ) = F( i, p ) * F( j, q );

  throwExceptionOnFailure( checkIfEqual( dTau_dPK2, expected_dTau_dPK2, 1e-12 ),
                           "FirstOrderDerived dTau_dPK2 with Diagonal Tensor mismatch in " +
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
                           "FirstOrderDerived dTau_dF with Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testKirchhoffStressFromPK2FirstOrderDerived_offDiagonalTensor()
{
  using namespace Fastor;

  double eps = 1e-8; // finite difference step

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
  auto [tau, dTau_dPK2, dTau_dF] = FirstOrderDerived::KirchhoffStressFromPK2( PK2, F );

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
                           "FirstOrderDerived tau with Off-Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // Test 2: Validate dTau/dPK2 analytically
  Tensor3333d ref_dTau_dPK2( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int p = 0; p < 3; ++p )
        for ( int q = 0; q < 3; ++q )
          ref_dTau_dPK2( i, j, p, q ) = F( i, p ) * F( j, q );

  throwExceptionOnFailure( checkIfEqual( dTau_dPK2, ref_dTau_dPK2, 1e-12 ),
                           "FirstOrderDerived dTau_dPK2 with Off-Diagonal Tensor mismatch in " +
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
                           "FirstOrderDerived dTau_dF with Off-Diagonal Tensor mismatch in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  // Register all test cases
  auto tests = std::vector< std::function< void() > >{ testKirchhoffStressFromPK2_identityTensor,
                                                       testKirchhoffStressFromPK2_diagonalTensor,
                                                       testKirchhoffStressFromPK2_offDiagonalTensor,
                                                       testKirchhoffStressFromPK2FirstOrderDerived_diagonalTensor,
                                                       testKirchhoffStressFromPK2FirstOrderDerived_offDiagonalTensor };

  // Execute and report any exceptions from failed tests
  executeTestsAndCollectExceptions( tests );
  return 0;
}
