#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Dense>
#include <cmath>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::NumericalAlgorithms;

// ---------------------------------------------------------------------------------------------
// No callers anywhere in the codebase currently use this class (HughesWingetWrapper in
// MarmotMaterialHughesWinget.h reimplements the same corotational algorithm independently rather
// than delegating to it), so there is no production usage pattern to reuse for these tests --
// verified instead against physical/algebraic identities and, for the tangent helpers, against
// the fixed projector tensors from MarmotKinematics.h that they are built from.
// ---------------------------------------------------------------------------------------------

void testPureStretchHasNoRotation()
{
  // FOld = I, FNew = a small symmetric (diagonal) stretch: no rotation should be detected.
  Eigen::Matrix3d FOld = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d FNew = Eigen::Matrix3d::Identity();
  FNew( 0, 0 )         = 1.002;
  FNew( 1, 1 )         = 0.999;

  HughesWinget hw( FOld, FNew, HughesWinget::AbaqusLike );

  throwExceptionOnFailure( checkIfEqual< double >( hw.getRotationIncrement(),
                                                   Eigen::Matrix3d::Identity().eval(),
                                                   1e-12 ),
                           "A pure stretch (no rotation) must give an identity rotation increment in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // The strain increment (Voigt) must match the documented formula (l = (FNew-FOld)*FMidStep^-1,
  // dEps = symmetric part of l, in Voigt notation), recomputed here independently.
  const Eigen::Matrix3d  FMidStep     = 0.5 * ( FNew + FOld );
  const Eigen::Matrix3d  lMat         = ( FNew - FOld ) * FMidStep.inverse();
  const Eigen::Matrix3d  dEpsMat      = 0.5 * ( lMat + lMat.transpose() );
  const Marmot::Vector6d dEpsExpected = ContinuumMechanics::VoigtNotation::voigtFromStrainMatrix< 3 >( dEpsMat );

  Marmot::Vector6d dEps = hw.getStrainIncrement();
  throwExceptionOnFailure( checkIfEqual< double >( dEps, dEpsExpected, 1e-12 ),
                           "Pure-stretch strain increment does not match the documented l = "
                           "(FNew-FOld)*FMidStep^-1 formula in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // With no rotation, rotateTensor() must be the identity map.
  Marmot::Vector6d someStress = Marmot::Vector6d::Zero();
  someStress << 10., 20., 30., 4., 5., 6.;
  throwExceptionOnFailure( checkIfEqual< double >( hw.rotateTensor( someStress ), someStress, 1e-10 ),
                           "rotateTensor() must be the identity map when there is no rotation increment in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testPureRotationPreservesIsotropicStress()
{
  // FOld = I, FNew = a small rotation about the z-axis: no stretching should be detected, and an
  // isotropic (spherical) stress -- invariant under any rotation -- must be returned unchanged.
  const double    phi  = 0.01; // rad
  Eigen::Matrix3d FOld = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d FNew = Eigen::Matrix3d::Identity();
  FNew( 0, 0 )         = std::cos( phi );
  FNew( 0, 1 )         = -std::sin( phi );
  FNew( 1, 0 )         = std::sin( phi );
  FNew( 1, 1 )         = std::cos( phi );

  HughesWinget hw( FOld, FNew, HughesWinget::AbaqusLike );

  throwExceptionOnFailure( checkIfEqual< double >( hw.getStrainIncrement(), Marmot::Vector6d::Zero().eval(), 1e-10 ),
                           "A pure rotation must give a zero strain increment in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // For a small rotation, the Hughes-Winget (Cayley-transform) incremental rotation must
  // approximate the exact rotation matrix.
  throwExceptionOnFailure( checkIfEqual< double >( hw.getRotationIncrement(), FNew, 1e-4 ),
                           "The incremental rotation does not approximate the applied small rotation in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  Marmot::Vector6d isotropicStress = Marmot::Vector6d::Zero();
  isotropicStress( 0 ) = isotropicStress( 1 ) = isotropicStress( 2 ) = -100.0; // hydrostatic pressure
  throwExceptionOnFailure( checkIfEqual< double >( hw.rotateTensor( isotropicStress ), isotropicStress, 1e-8 ),
                           "An isotropic stress state must be invariant under rotateTensor() in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testRotateTensorMatchesDirectSimilarityTransform()
{
  // General (non-trivial) FOld/FNew, so dR is neither the identity nor a simple z-rotation.
  Eigen::Matrix3d FOld = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d FNew( 3, 3 );
  // clang-format off
  FNew << 1.02,  0.03, -0.01,
          0.04,  0.98,  0.02,
         -0.02,  0.01,  1.01;
  // clang-format on

  HughesWinget hw( FOld, FNew, HughesWinget::AbaqusLike );

  Marmot::Vector6d stress = Marmot::Vector6d::Zero();
  stress << 50., -20., 10., 5., -3., 2.;

  const Marmot::Vector6d rotated = hw.rotateTensor( stress );

  // Recompute the same similarity transform independently via the plain 3x3 matrices, rather than
  // by re-deriving rotateTensor()'s own Voigt-conversion internals.
  const Eigen::Matrix3d  dR                 = hw.getRotationIncrement();
  const Eigen::Matrix3d  stressMat          = ContinuumMechanics::VoigtNotation::voigtToStress< double >( stress );
  const Eigen::Matrix3d  rotatedMatExpected = dR * stressMat * dR.transpose();
  const Marmot::Vector6d rotatedExpected    = ContinuumMechanics::VoigtNotation::stressToVoigt< double >(
    rotatedMatExpected );

  throwExceptionOnFailure( checkIfEqual< double >( rotated, rotatedExpected, 1e-10 ),
                           "rotateTensor() does not match dR * stress * dR^T in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // dR must be a proper rotation (orthogonal, unit determinant), as returned by a Cayley transform
  // of a skew matrix.
  throwExceptionOnFailure( checkIfEqual< double >( ( dR.transpose() * dR ).eval(),
                                                   Eigen::Matrix3d::Identity().eval(),
                                                   1e-10 ) &&
                             checkIfEqual( dR.determinant(), 1.0, 1e-10 ),
                           "The incremental rotation matrix must be orthogonal with unit determinant in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// compute_dS_dF() / compute_dScalar_dF(): both are, for fixed inputs, purely linear combinations
// of the caller-supplied tensors with the *fixed* (state-independent) projector tensors
// dOmega_dVelocityGradient / dStretchingRate_dVelocityGradient from MarmotKinematics.h. Choosing
// FOld = FNew = I makes FInv = I and isolates each term in turn (zero stress isolates the
// Jaumann/material term; zero dChauchydEps isolates the rotational term), so both can be checked
// against a closed-form reference built directly from those same public projector tensors --
// without needing to re-derive or duplicate the (approximate, Abaqus-UMAT-style) tangent formula.
// ---------------------------------------------------------------------------------------------

void testComputeDScalarDFMatchesProjectorContraction()
{
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  HughesWinget    hw( I, I, HughesWinget::AbaqusLike );

  Marmot::Vector6d dScalarDEps = Marmot::Vector6d::Zero();
  dScalarDEps << 1., 2., 3., 4., 5., 6.;

  const Eigen::Matrix3d result = hw.compute_dScalar_dF( I, dScalarDEps );

  using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;
  Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ )
      for ( int ij = 0; ij < 6; ij++ )
        expected( k, l ) += dScalarDEps( ij ) * dStretchingRate_dVelocityGradient( ij, k, l );

  throwExceptionOnFailure( checkIfEqual< double >( result, expected, 1e-12 ),
                           "compute_dScalar_dF() does not match the direct contraction with "
                           "dStretchingRate_dVelocityGradient in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputeDSDFJaumannTermMatchesProjectorContraction()
{
  // Zero stress isolates the Jaumann (material-tangent-driven) term.
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  HughesWinget    hw( I, I, HughesWinget::AbaqusLike );

  Marmot::Vector6d zeroStress = Marmot::Vector6d::Zero();

  Marmot::Matrix6d dChauchydEps = Marmot::Matrix6d::Zero();
  for ( int i = 0; i < 6; i++ )
    for ( int j = 0; j < 6; j++ )
      dChauchydEps( i, j ) = 1.0 + i + 2.0 * j;

  const EigenTensors::Tensor633d result = hw.compute_dS_dF( zeroStress, I, dChauchydEps );

  using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;
  bool matches = true;
  for ( int ij = 0; ij < 6 && matches; ij++ )
    for ( int k = 0; k < 3 && matches; k++ )
      for ( int l = 0; l < 3 && matches; l++ ) {
        double expected = 0.0;
        for ( int mn = 0; mn < 6; mn++ )
          expected += dChauchydEps( ij, mn ) * dStretchingRate_dVelocityGradient( mn, k, l );
        matches = matches && checkIfEqual( result( ij, k, l ), expected, 1e-10 );
      }

  throwExceptionOnFailure( matches,
                           "compute_dS_dF() (Jaumann term, zero stress) does not match the direct contraction "
                           "with dStretchingRate_dVelocityGradient in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testComputeDSDFRotationalTermMatchesProjectorContraction()
{
  // Zero dChauchydEps isolates the rotational (stress-spin) term.
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  HughesWinget    hw( I, I, HughesWinget::AbaqusLike );

  Marmot::Vector6d stress = Marmot::Vector6d::Zero();
  stress << 10., 20., 30., 4., 5., 6.;
  const Eigen::Matrix3d stressMat = ContinuumMechanics::VoigtNotation::voigtToStress< double >( stress );

  Marmot::Matrix6d zeroTangent = Marmot::Matrix6d::Zero();

  const EigenTensors::Tensor633d result = hw.compute_dS_dF( stress, I, zeroTangent );

  using namespace Marmot::ContinuumMechanics::Kinematics::VelocityGradient;
  using namespace Marmot::ContinuumMechanics::TensorUtility;
  bool matches = true;
  for ( int ij = 0; ij < 6 && matches; ij++ ) {
    auto [i, j] = IndexNotation::fromVoigt< 3 >( ij );
    for ( int k = 0; k < 3 && matches; k++ )
      for ( int l = 0; l < 3 && matches; l++ ) {
        double expected = 0.0;
        for ( int m = 0; m < 3; m++ )
          expected += dOmega_dVelocityGradient( i, m, k, l ) * stressMat( m, j ) +
                      dOmega_dVelocityGradient( j, m, k, l ) * stressMat( i, m );
        matches = matches && checkIfEqual( result( ij, k, l ), expected, 1e-10 );
      }
  }

  throwExceptionOnFailure( matches,
                           "compute_dS_dF() (rotational term, zero tangent) does not match the direct "
                           "contraction with dOmega_dVelocityGradient in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testPureStretchHasNoRotation,
    testPureRotationPreservesIsotropicStress,
    testRotateTensorMatchesDirectSimilarityTransform,
    testComputeDScalarDFMatchesProjectorContraction,
    testComputeDSDFJaumannTermMatchesProjectorContraction,
    testComputeDSDFRotationalTermMatchesProjectorContraction,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
