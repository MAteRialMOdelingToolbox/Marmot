#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics;
using namespace Marmot::ContinuumMechanics::DeformationMeasures;

// Expected C = F^T * F
Tensor33d C_expected( const FastorStandardTensors::Tensor33d& F_ )
{
  FastorStandardTensors::Tensor33d temp( 0.0 );
  for ( int I = 0; I < 3; ++I )
    for ( int J = 0; J < 3; ++J )
      for ( int K = 0; K < 3; ++K )
        temp( I, J ) += F_( K, I ) * F_( K, J );
  return temp;
};

// Expected derivative dC_IJ / dF_KL = delta_IL * F_KJ + F_KI * delta_JL
Tensor3333d dC_dF_expected( const FastorStandardTensors::Tensor33d& F_ )
{
  FastorStandardTensors::Tensor3333d temp( 0.0 );
  for ( int I = 0; I < 3; ++I )
    for ( int J = 0; J < 3; ++J )
      for ( int K = 0; K < 3; ++K )
        for ( int L = 0; L < 3; ++L ) {
          int delta_IL       = ContinuumMechanics::TensorUtility::d( I, L );
          int delta_JL       = ContinuumMechanics::TensorUtility::d( J, L );
          temp( I, J, K, L ) = delta_IL * F_( K, J ) + F_( K, I ) * delta_JL;
        }
  return temp;
};

void test_rightCauchyGreen()
{
  // Testcase 1: Identity matrix: F = I  -->  C = I
  {
    FastorStandardTensors::Tensor33t< double > F;
    // clang-format off
    F = { {1.,0.,0.},
          {0.,1.,0.},
          {0.,0.,1.} };
    // clang-format on

    auto C_computed = DeformationMeasures::rightCauchyGreen( F );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
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

    auto C_computed = DeformationMeasures::rightCauchyGreen( F );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
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

    auto C_computed = DeformationMeasures::rightCauchyGreen( F );

    throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3" );
  }
}

void test_rightCauchyGreenirstOrderDerived()
{
  // Testcase 1: Identity matrix F
  FastorStandardTensors::Tensor33t< double > F;
  F = { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } };

  auto [C_computed, dC_dF_computed] = FirstOrderDerived::rightCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1 (computation of C)" );

  throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 1 (computation of dC_dF)" );

  // Testcase 2: Diagonal matrix F
  F = { { 2., 0., 0. }, { 0., 3., 0. }, { 0., 0., 4. } };

  std::tie( C_computed, dC_dF_computed ) = FirstOrderDerived::rightCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2 (computation of C)" );

  throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 2 (computation of dC_dF)" );

  // Testcase 3: Non-diagonal matrix F
  F = { { 1., 1., 1. }, { 0., 1., 1. }, { 0., 0., 1. } };

  std::tie( C_computed, dC_dF_computed ) = FirstOrderDerived::rightCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( C_computed, C_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3 (computation of C)" );

  throwExceptionOnFailure( checkIfEqual( dC_dF_computed, dC_dF_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for Testcase 3 (computation of dC_dF)" );
}

Tensor33d b_expected( const FastorStandardTensors::Tensor33d& F_ )
{
  FastorStandardTensors::Tensor33d temp( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int K = 0; K < 3; ++K )
        temp( i, j ) += F_( i, K ) * F_( j, K );
  return temp;
};

Tensor3333d db_dF_expected( const FastorStandardTensors::Tensor33d& F_ )
{
  FastorStandardTensors::Tensor3333d temp( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int K = 0; K < 3; ++K ) {
          int delta_ik       = ContinuumMechanics::TensorUtility::d( i, k );
          int delta_jk       = ContinuumMechanics::TensorUtility::d( j, k );
          temp( i, j, k, K ) = delta_ik * F_( j, K ) + F_( i, K ) * delta_jk;
        }
  return temp;
};

void test_leftCauchyGreen()
{
  // Simple test: for F = I, B = I
  FastorStandardTensors::Tensor33d F = FastorStandardTensors::Spatial3D::I;
  auto                             b = leftCauchyGreen( F );
  throwExceptionOnFailure( checkIfEqual( b, F ), "test_leftCauchGreen failed for zero deformation" );

  // Simple test: for F = 2*I, B = 4*I
  F *= 2.0;
  b = leftCauchyGreen( F );
  throwExceptionOnFailure( checkIfEqual( b, Tensor33d( 2. * F ) ), "failed for uniform expansion" );

  // Simple test: for F = diag(1,2,3), B = diag(1,4,9)
  F         = FastorStandardTensors::Spatial3D::I;
  F( 1, 1 ) = 2.0;
  F( 2, 2 ) = 3.0;
  b         = leftCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( b, b_expected( F ) ), "failed for non-uniform expansion" );

  // test for simple shear
  F         = FastorStandardTensors::Spatial3D::I;
  F( 0, 1 ) = 1.0; // simple shear in x-y plane
  b         = leftCauchyGreen( F );
  throwExceptionOnFailure( checkIfEqual( b, b_expected( F ) ), "failed for simple shear" );
}

void test_leftCauchyGreenFirstOrderDerived()
{
  // Testcase 1: Identity matrix F
  FastorStandardTensors::Tensor33t< double > F;
  F = { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } };

  auto [b_computed, db_dF_computed] = FirstOrderDerived::leftCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( b_computed, b_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for zero deformation (computation of b)" );

  // Expected derivative db_ij / dF_kK = delta_ik * F_jK + F_iK * delta_jk

  throwExceptionOnFailure( checkIfEqual( db_dF_computed, db_dF_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " failed for zero deformation (computation of db_dF)" );

  // simple shear test
  F                                             = FastorStandardTensors::Spatial3D::I;
  F( 0, 1 )                                     = 1.0; // simple shear in x-y
  auto [b_shear_computed, db_dF_shear_computed] = FirstOrderDerived::leftCauchyGreen( F );

  throwExceptionOnFailure( checkIfEqual( b_shear_computed, b_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for simple shear (computation of b)" );

  throwExceptionOnFailure( checkIfEqual( db_dF_shear_computed, db_dF_expected( F ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for simple shear (computation of db_dF)" );
}
int main()
{
  auto tests = std::vector< std::function< void() > >{ test_rightCauchyGreen,
                                                       test_rightCauchyGreenirstOrderDerived,
                                                       test_leftCauchyGreen };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
