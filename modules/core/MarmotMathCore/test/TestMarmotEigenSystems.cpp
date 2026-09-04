#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"

using Marmot::MakeString;
using namespace Marmot::Testing;
using namespace Marmot::TensorUtility::FastorTensors::StandardTensors;
using namespace Marmot::Math;
using namespace Marmot::NumericalAlgorithms::AutomaticDifferentiation;

// Helper: check that V^T * V ≈ I, i.e. the eigenvector columns are orthonormal
static void checkOrthonormality( const Tensor33d& V, const std::string& caller, double tol = 1e-10 )
{
  for ( int k = 0; k < 3; k++ ) {
    for ( int l = 0; l < 3; l++ ) {
      double vtv_kl = 0.0;
      for ( int i = 0; i < 3; i++ )
        vtv_kl += V( i, k ) * V( i, l );
      const double expected = ( k == l ) ? 1.0 : 0.0;
      throwExceptionOnFailure( checkIfEqual( vtv_kl, expected, tol ),
                               MakeString() << caller << " eigenvectors not orthonormal at (" << k << "," << l
                                            << "): got " << vtv_kl << " expected " << expected );
    }
  }
}

// Helper: check spectral reconstruction A ≈ V * diag(eigenvalues) * V^T
static void checkReconstruction( const Tensor33d&   A,
                                 const Tensor3d&    eigenvalues,
                                 const Tensor33d&   V,
                                 const std::string& caller,
                                 double             tol = 1e-10 )
{
  Tensor33d A_reconstructed( 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        A_reconstructed( i, j ) += V( i, k ) * eigenvalues( k ) * V( j, k );

  throwExceptionOnFailure( checkIfEqual< double >( A_reconstructed, A, tol ),
                           MakeString() << caller << " spectral reconstruction failed" );
}

// ============================================================
//  Tests for computeEigenSystemJacobi
// ============================================================

void testJacobi_identity()
{
  Tensor33d A;
  A.eye();

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  for ( int k = 0; k < 3; k++ )
    throwExceptionOnFailure( checkIfEqual( eigenvalues( k ), 1.0, 1e-12 ),
                             MakeString() << __PRETTY_FUNCTION__ << " eigenvalue " << k << " != 1" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testJacobi_diagonal()
{
  // Diagonal matrix with distinct eigenvalues already in descending order
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 5.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 1.0;

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 5.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 3.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 1.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testJacobi_symmetric_distinct_eigenvalues()
{
  // A = [[5, 1, 0], [1, 5, 0], [0, 0, 3]]
  // 2x2 block eigenvalues: (5-λ)^2 - 1 = 0 → λ = 4 or 6
  // Sorted descending eigenvalues: [6, 4, 3]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 5.0;
  A( 0, 1 ) = 1.0;
  A( 1, 0 ) = 1.0;
  A( 1, 1 ) = 5.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 6.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 4.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 3.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testJacobi_fully_dense_symmetric()
{
  // A = [[3, 1, 1], [1, 3, 1], [1, 1, 3]]
  // Characteristic polynomial: (λ-5)(λ-2)^2 = 0 → eigenvalues [5, 2, 2]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 3.0;
  A( 0, 1 ) = 1.0;
  A( 0, 2 ) = 1.0;
  A( 1, 0 ) = 1.0;
  A( 1, 1 ) = 3.0;
  A( 1, 2 ) = 1.0;
  A( 2, 0 ) = 1.0;
  A( 2, 1 ) = 1.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 5.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 2.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 2.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testJacobi_negative_eigenvalues()
{
  // Diagonal matrix with negative entries, provided in non-descending order
  // diag(-1, -2, 3) → sorted descending: [3, -1, -2]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = -1.0;
  A( 1, 1 ) = -2.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 3.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), -1.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), -2.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testJacobi_sorted_descending()
{
  // Provide a diagonal matrix with eigenvalues in non-descending order
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 1.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 2.0;

  auto [eigenvalues, V] = computeEigenSystemJacobi( A );

  throwExceptionOnFailure( eigenvalues( 0 ) >= eigenvalues( 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << " eigenvalues not sorted descending (λ₀ >= λ₁)" );
  throwExceptionOnFailure( eigenvalues( 1 ) >= eigenvalues( 2 ),
                           MakeString() << __PRETTY_FUNCTION__ << " eigenvalues not sorted descending (λ₁ >= λ₂)" );
}

void testJacobi_automatic_differentiation()
{
  // Diagonal matrix diag(3, 2, 1). For this matrix λ_k = A(k,k), so
  // dλ_0/dA has only one non-zero entry: dλ_0/dA(0,0) = 1
  Tensor33d A_double( 0.0 );
  A_double( 0, 0 ) = 3.0;
  A_double( 1, 1 ) = 2.0;
  A_double( 2, 2 ) = 1.0;

  std::function< dual( const Tensor33t< dual >& ) > f0 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemJacobi< dual >( A_d ).first( 0 );
  };
  Tensor33d dlambda0_dA = df_dT( f0, A_double );

  Tensor33d expected0( 0.0 );
  expected0( 0, 0 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda0_dA, expected0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₀ w.r.t. A failed" );

  std::function< dual( const Tensor33t< dual >& ) > f1 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemJacobi< dual >( A_d ).first( 1 );
  };
  Tensor33d dlambda1_dA = df_dT( f1, A_double );

  Tensor33d expected1( 0.0 );
  expected1( 1, 1 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda1_dA, expected1, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₁ w.r.t. A failed" );

  std::function< dual( const Tensor33t< dual >& ) > f2 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemJacobi< dual >( A_d ).first( 2 );
  };
  Tensor33d dlambda2_dA = df_dT( f2, A_double );

  Tensor33d expected2( 0.0 );
  expected2( 2, 2 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda2_dA, expected2, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₂ w.r.t. A failed" );
}

// ============================================================
//  Tests for computeEigenSystemWithCardano
// ============================================================

void testCardano_identity()
{
  Tensor33d A;
  A.eye();

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  for ( int k = 0; k < 3; k++ )
    throwExceptionOnFailure( checkIfEqual( eigenvalues( k ), 1.0, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " eigenvalue " << k << " != 1" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_diagonal()
{
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 5.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 1.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 5.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 3.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 1.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_symmetric_distinct_eigenvalues()
{
  // A = [[5, 1, 0], [1, 5, 0], [0, 0, 3]], sorted eigenvalues: [6, 4, 3]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 5.0;
  A( 0, 1 ) = 1.0;
  A( 1, 0 ) = 1.0;
  A( 1, 1 ) = 5.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 6.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 4.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 3.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_fully_dense_symmetric()
{
  // A = [[3, 1, 1], [1, 3, 1], [1, 1, 3]], eigenvalues [5, 2, 2]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 3.0;
  A( 0, 1 ) = 1.0;
  A( 0, 2 ) = 1.0;
  A( 1, 0 ) = 1.0;
  A( 1, 1 ) = 3.0;
  A( 1, 2 ) = 1.0;
  A( 2, 0 ) = 1.0;
  A( 2, 1 ) = 1.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 5.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 2.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 2.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_all_repeated()
{
  // Isotropic case: all eigenvalues equal. Cardano has a dedicated branch for this.
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 4.0;
  A( 1, 1 ) = 4.0;
  A( 2, 2 ) = 4.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  for ( int k = 0; k < 3; k++ )
    throwExceptionOnFailure( checkIfEqual( eigenvalues( k ), 4.0, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " eigenvalue " << k << " wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_two_repeated_top()
{
  // λ₀ = λ₁ = 3, λ₂ = 1 (root_01_repeated branch in Cardano)
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 3.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 1.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  // Cardano uses trigonometric functions; slightly looser tolerance for repeated roots
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 3.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 3.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 1.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__, 1e-8 );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__, 1e-8 );
}

void testCardano_two_repeated_bottom()
{
  // λ₀ = 5, λ₁ = λ₂ = 2 (root_12_repeated branch in Cardano)
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 5.0;
  A( 1, 1 ) = 2.0;
  A( 2, 2 ) = 2.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  // Cardano uses trigonometric functions; slightly looser tolerance for repeated roots
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 5.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), 2.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), 2.0, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__, 1e-8 );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__, 1e-8 );
}

void testCardano_negative_eigenvalues()
{
  // diag(-1, -2, 3) → sorted descending: [3, -1, -2]
  Tensor33d A( 0.0 );
  A( 0, 0 ) = -1.0;
  A( 1, 1 ) = -2.0;
  A( 2, 2 ) = 3.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  throwExceptionOnFailure( checkIfEqual( eigenvalues( 0 ), 3.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " first eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 1 ), -1.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " second eigenvalue wrong" );
  throwExceptionOnFailure( checkIfEqual( eigenvalues( 2 ), -2.0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " third eigenvalue wrong" );

  checkOrthonormality( V, __PRETTY_FUNCTION__ );
  checkReconstruction( A, eigenvalues, V, __PRETTY_FUNCTION__ );
}

void testCardano_sorted_descending()
{
  // Provide matrix with eigenvalues in non-descending order on diagonal
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 1.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 2.0;

  auto [eigenvalues, V] = computeEigenSystemWithCardano( A );

  throwExceptionOnFailure( eigenvalues( 0 ) >= eigenvalues( 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << " eigenvalues not sorted descending (λ₀ >= λ₁)" );
  throwExceptionOnFailure( eigenvalues( 1 ) >= eigenvalues( 2 ),
                           MakeString() << __PRETTY_FUNCTION__ << " eigenvalues not sorted descending (λ₁ >= λ₂)" );
}

void testCardano_automatic_differentiation()
{
  // Diagonal matrix diag(3, 2, 1). For this matrix λ_k = A(k,k), so
  // dλ_0/dA has only one non-zero entry: dλ_0/dA(0,0) = 1
  Tensor33d A_double( 0.0 );
  A_double( 0, 0 ) = 3.0;
  A_double( 1, 1 ) = 2.0;
  A_double( 2, 2 ) = 1.0;

  std::function< dual( const Tensor33t< dual >& ) > f0 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemWithCardano< dual >( A_d ).first( 0 );
  };
  Tensor33d dlambda0_dA = df_dT( f0, A_double );

  Tensor33d expected0( 0.0 );
  expected0( 0, 0 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda0_dA, expected0, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₀ w.r.t. A failed" );

  std::function< dual( const Tensor33t< dual >& ) > f2 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemWithCardano< dual >( A_d ).first( 2 );
  };
  Tensor33d dlambda2_dA = df_dT( f2, A_double );

  Tensor33d expected2( 0.0 );
  expected2( 2, 2 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda2_dA, expected2, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₂ w.r.t. A failed" );

  std::function< dual( const Tensor33t< dual >& ) > f1 = []( const Tensor33t< dual >& A_d ) {
    return computeEigenSystemWithCardano< dual >( A_d ).first( 1 );
  };
  Tensor33d dlambda1_dA = df_dT( f1, A_double );

  Tensor33d expected1( 0.0 );
  expected1( 1, 1 ) = 1.0;
  throwExceptionOnFailure( checkIfEqual< double >( dlambda1_dA, expected1, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " gradient of λ₁ w.r.t. A failed" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testJacobi_identity,
    testJacobi_diagonal,
    testJacobi_symmetric_distinct_eigenvalues,
    testJacobi_fully_dense_symmetric,
    testJacobi_negative_eigenvalues,
    testJacobi_sorted_descending,
    testJacobi_automatic_differentiation,
    testCardano_identity,
    testCardano_diagonal,
    testCardano_symmetric_distinct_eigenvalues,
    testCardano_fully_dense_symmetric,
    testCardano_all_repeated,
    testCardano_two_repeated_top,
    testCardano_two_repeated_bottom,
    testCardano_negative_eigenvalues,
    testCardano_sorted_descending,
    testCardano_automatic_differentiation,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
