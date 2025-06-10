#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics;
using namespace Marmot::ContinuumMechanics::DeformationMeasures;

void test_CauchyGreen()
{
  // Testcase 1: Identity matrix: F = I  -->  C = I
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
    F = { {1.,0.,0.},
          {0.,1.,0.},
          {0.,0.,1.} };
    // clang-format on

    auto C_computed = DeformationMeasures::CauchyGreen( F );

    FastorStandardTensors::Tensor33t< double > C_expected;
    // clang-format off
    C_expected = { {1.,0.,0.},
                   {0.,1.,0.},
                   {0.,0.,1.} };
    // clang-format on

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1" );
  }

  // Testcase 2: Diagonal matrix: F = diag(2,3,4)  -->  C = diag(4,9,16)
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
    F = { {2.,0.,0.},
          {0.,3.,0.},
          {0.,0.,4.} };
    // clang-format on

    auto C_computed = DeformationMeasures::CauchyGreen( F );

    FastorStandardTensors::Tensor33t< double > C_expected;
    // clang-format off
    C_expected = { {4.,0.,0.},
                   {0.,9.,0.},
                   {0.,0.,16.} };
    // clang-format on

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2" );
  }

  // Testcase 3: Non-diagonal matrix:
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
    F = { {1., 1., 1.},
          {0., 1., 1.},
          {0., 0., 1.} };
    // clang-format on

    auto C_computed = DeformationMeasures::CauchyGreen( F );

    // Expected: C = F^T * F = [ [1,0,0; 1,1,0; 1,1,1] * [1,1,1; 0,1,1; 0,0,1] ]
    FastorStandardTensors::Tensor33t< double > C_expected;
    // clang-format off
    C_expected = { {1., 1., 1.},
                   {1., 2., 2.},
                   {1., 2., 3.} };
    // clang-format on

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }
}

void test_FirstOrderDerived()
{
  // Testcase 1: Identity matrix F
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
		F = { {1.,0.,0.},
          {0.,1.,0.},
          {0.,0.,1.} };
    // clang-format on

    auto [C_computed, dC_dF_computed] = FirstOrderDerived::CauchyGreen( F );

    // C = F^T * F -> when F is diagonal: C_II = F_II^2
    FastorStandardTensors::Tensor33t< double > C_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      C_expected( I, I ) = F( 0, I ) * F( 0, I ) + F( 1, I ) * F( 1, I ) + F( 2, I ) * F( 2, I );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1 (computation of C)" );

    // Expected derivative dC_IJ / dF_KL = delta_IL * F_KJ + F_KI * delta_JL
    FastorStandardTensors::Tensor3333t< double > dC_dF_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      for ( int J = 0; J < 3; ++J )
        for ( int K = 0; K < 3; ++K )
          for ( int L = 0; L < 3; ++L ) {
            double delta_IL              = ( I == L ? 1.0 : 0.0 );
            double delta_JL              = ( J == L ? 1.0 : 0.0 );
            dC_dF_expected( I, J, K, L ) = delta_IL * F( K, J ) + F( K, I ) * delta_JL;
          }

    throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1 (computation of dC_dF)" );
  }

  // Testcase 2: Diagonal matrix F
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
		F = { {2.,0.,0.},
		      {0.,3.,0.},
		      {0.,0.,4.} };
    // clang-format on

    auto [C_computed, dC_dF_computed] = FirstOrderDerived::CauchyGreen( F );

    // C = F^T * F -> when F is diagonal: C_II = F_II^2
    FastorStandardTensors::Tensor33t< double > C_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      C_expected( I, I ) = F( 0, I ) * F( 0, I ) + F( 1, I ) * F( 1, I ) + F( 2, I ) * F( 2, I );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2 (computation of C)" );

    // Expected derivative dC_IJ / dF_KL = delta_IL * F_KJ + F_KI * delta_JL
    FastorStandardTensors::Tensor3333t< double > dC_dF_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      for ( int J = 0; J < 3; ++J )
        for ( int K = 0; K < 3; ++K )
          for ( int L = 0; L < 3; ++L ) {
            double delta_IL              = ( I == L ? 1.0 : 0.0 );
            double delta_JL              = ( J == L ? 1.0 : 0.0 );
            dC_dF_expected( I, J, K, L ) = delta_IL * F( K, J ) + F( K, I ) * delta_JL;
          }

    throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2 (computation of dC_dF)" );
  }

  // Testcase 3: Non-diagonal matrix F
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
		F = { {1.,1.,1.},
		      {0.,1.,1.},
		      {0.,0.,1.} };
    // clang-format on

    auto [C_computed, dC_dF_computed] = FirstOrderDerived::CauchyGreen( F );

    // C = F^T * F = sum_i F_iI * F_iJ
    FastorStandardTensors::Tensor33t< double > C_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      for ( int J = 0; J < 3; ++J )
        for ( int i = 0; i < 3; ++i )
          C_expected( I, J ) += F( i, I ) * F( i, J );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3 (computation of C)" );

    // dC_IJ / dF_KL = delta_IL * F_KJ + F_KI * delta_JL
    FastorStandardTensors::Tensor3333t< double > dC_dF_expected( 0.0 );
    for ( int I = 0; I < 3; ++I )
      for ( int J = 0; J < 3; ++J )
        for ( int K = 0; K < 3; ++K )
          for ( int L = 0; L < 3; ++L ) {
            double delta_IL              = ( I == L ? 1.0 : 0.0 );
            double delta_JL              = ( J == L ? 1.0 : 0.0 );
            dC_dF_expected( I, J, K, L ) = delta_IL * F( K, J ) + F( K, I ) * delta_JL;
          }

    throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3 (computation of dC_dF)" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ test_CauchyGreen, test_FirstOrderDerived };

  executeTestsAndCollectExceptions( tests );
  return 0;
}