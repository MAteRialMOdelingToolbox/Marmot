#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <algorithm>
#include <cmath>

using namespace Marmot::Testing;

using namespace Marmot;

void testStrainToVoigt()
{
  Eigen::Matrix< double, 3, 3 > strain;
  strain << 1, 4, 5, 4, 2, 6, 5, 6, 3;
  Vector6d   strainVoigtGold = { 1, 2, 3, 4 * 2, 5 * 2, 6 * 2 };
  const auto strainVoigt     = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt( strain );
  throwExceptionOnFailure( checkIfEqual< double >( strainVoigt, strainVoigtGold ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testVoigtToPlaneVoigt()
{
  Vector6d   strainVoigt = { 1, 2, 3, 4, 5, 6 };
  const auto planeVoigt  = Marmot::ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( strainVoigt );
  throwExceptionOnFailure( checkIfEqual< double >( planeVoigt, Vector3d( 1, 2, 4 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testPrincipalStrains()
{
  Vector6d   strainVoigt      = { 1, 1, 1, 0, 0, 0 };
  const auto principalStrains = Marmot::ContinuumMechanics::VoigtNotation::Invariants::principalStrains( strainVoigt );
  throwExceptionOnFailure( checkIfEqual< double >( principalStrains, Vector3d( 1, 1, 1 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testPrincipalStresses()
{
  Vector6d   stressVoigt       = { 1, 1, 1, 0, 0, 0 };
  const auto principalStresses = Marmot::ContinuumMechanics::VoigtNotation::Invariants::principalStresses(
    stressVoigt );
  throwExceptionOnFailure( checkIfEqual< double >( principalStresses, Vector3d( 1, 1, 1 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testSortedPrincipalStrains()
{
  Vector6d   strainVoigt      = { 1, 3, 2, 0, 0, 0 };
  const auto principalStrains = Marmot::ContinuumMechanics::VoigtNotation::Invariants::sortedPrincipalStrains(
    strainVoigt );
  throwExceptionOnFailure( checkIfEqual< double >( principalStrains, Vector3d( 3, 2, 1 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testPrincipalStressDirections()
{

  Vector6d stressVoigt = { 1, 1, 1, 0, 0, 0 };
  const auto
    principalStressDirections = Marmot::ContinuumMechanics::VoigtNotation::Invariants::principalStressesDirections(
      stressVoigt );
  throwExceptionOnFailure( checkIfEqual( principalStressDirections.norm(), std::sqrt( 3. ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testVonMisesEquivalentStress()
{
  {
    Vector6d stressVoigt = { 1, 1, 1, 0, 0, 0 };
    const auto
      vonMisesEquivalentStress = Marmot::ContinuumMechanics::VoigtNotation::Invariants::vonMisesEquivalentStress(
        stressVoigt );
    throwExceptionOnFailure( checkIfEqual( vonMisesEquivalentStress, 0.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed" );
  }
  {
    Vector6d stressVoigt = { 1, 0, 0, 0, 0, 0 };
    const auto
      vonMisesEquivalentStress = Marmot::ContinuumMechanics::VoigtNotation::Invariants::vonMisesEquivalentStress(
        stressVoigt );
    throwExceptionOnFailure( checkIfEqual( vonMisesEquivalentStress, 1.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed" );
  }
}

void testVonMisesEquivalentStrain()
{
  {
    Vector6d strainVoigt = { 1, 1, 1, 0, 0, 0 };
    const auto
      vonMisesEquivalentStrain = Marmot::ContinuumMechanics::VoigtNotation::Invariants::vonMisesEquivalentStrain(
        strainVoigt );
    throwExceptionOnFailure( checkIfEqual( vonMisesEquivalentStrain, std::sqrt( 2.0 ) ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed" );
  }
  {
    Vector6d strainVoigt = { 1, 0, 0, 0, 0, 0 };
    const auto
      vonMisesEquivalentStrain = Marmot::ContinuumMechanics::VoigtNotation::Invariants::vonMisesEquivalentStrain(
        strainVoigt );
    throwExceptionOnFailure( checkIfEqual( vonMisesEquivalentStrain, std::sqrt( 2. / 3.0 ) ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed" );
  }
}

void testI1()
{
  const Vector6d stress = { 1, 2, 3, 0, 0, 0 };
  const auto     I1     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress );
  throwExceptionOnFailure( checkIfEqual( I1, 6 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testI2()
{
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };
  const auto     I2     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::I2( stress );
  throwExceptionOnFailure( checkIfEqual( I2, -66 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testI2Strain()
{
  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };
  const auto     I2     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::I2Strain( strain );
  throwExceptionOnFailure( checkIfEqual( I2, -8.25 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testI3()
{
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };
  const auto     I3     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::I3( stress );
  throwExceptionOnFailure( checkIfEqual( I3, ContinuumMechanics::VoigtNotation::voigtToStress( stress ).determinant() ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testI3Strain()
{
  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };
  const auto     I3     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::I3Strain( strain );
  throwExceptionOnFailure( checkIfEqual( I3, ContinuumMechanics::VoigtNotation::voigtToStrain( strain ).determinant() ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testJ2()
{
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };
  const auto     J2     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::J2( stress );
  throwExceptionOnFailure( checkIfEqual( J2,
                                         1.0 / 3.0 *
                                           ( Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) *
                                               Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) -
                                             3 *
                                               Marmot::ContinuumMechanics::VoigtNotation::Invariants::I2( stress ) ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testJ3()
{
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };
  const auto     J3     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::J3( stress );
  throwExceptionOnFailure( checkIfEqual( J3,
                                         2.0 / 27.0 *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) -
                                           1. / 3 *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( stress ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I2( stress ) +
                                           1 * Marmot::ContinuumMechanics::VoigtNotation::Invariants::I3( stress ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testJ3Strain()
{

  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };
  const auto     J3     = Marmot::ContinuumMechanics::VoigtNotation::Invariants::J3Strain( strain );
  throwExceptionOnFailure( checkIfEqual( J3,
                                         2.0 / 27.0 *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( strain ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( strain ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( strain ) -
                                           1. / 3 *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I1( strain ) *
                                             Marmot::ContinuumMechanics::VoigtNotation::Invariants::I2Strain( strain ) +
                                           1 * Marmot::ContinuumMechanics::VoigtNotation::Invariants::I3Strain(
                                                 strain ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dStressMean_dStress()
{

  using namespace Marmot::ContinuumMechanics::VoigtNotation::Derivatives;

  const Vector6d dStressMean_dStressGold = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0, 0, 0 };

  throwExceptionOnFailure( checkIfEqual( dStressMean_dStress().norm(), dStressMean_dStressGold.norm() ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dRho_dStress()
{

  using namespace Marmot::ContinuumMechanics::VoigtNotation::Derivatives;
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };

  const auto hw = ContinuumMechanics::HaighWestergaard::haighWestergaard( stress );

  const auto rho = []( const Vector6d& stress ) {
    Eigen::MatrixXd rho( 1, 1 );
    rho << ContinuumMechanics::HaighWestergaard::haighWestergaard( stress ).rho;
    return rho;
  };

  const auto dRho_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( rho, stress );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dRho_dStress( hw.rho,
                                                                                                               stress )
                                           .norm(),
                                         dRho_dStress_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dRhoStrain_dStrain()
{

  using namespace Marmot::ContinuumMechanics::VoigtNotation::Derivatives;
  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };

  const auto hw = ContinuumMechanics::HaighWestergaard::haighWestergaardFromStrain( strain );

  const auto rho = []( const Vector6d& strain ) {
    Eigen::MatrixXd rho( 1, 1 );
    rho << ContinuumMechanics::HaighWestergaard::haighWestergaardFromStrain( strain ).rho;
    return rho;
  };

  const auto dRho_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( rho, strain );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::
                                           dRhoStrain_dStrain( hw.rho, strain )
                                             .norm(),
                                         dRho_dStress_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dTheta_dStress()
{

  using namespace Marmot::ContinuumMechanics::VoigtNotation::Derivatives;
  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };

  const auto theta = []( const Vector6d& stress ) {
    Eigen::MatrixXd theta( 1, 1 );
    theta << ContinuumMechanics::HaighWestergaard::haighWestergaard( stress ).theta;
    return theta;
  };

  const auto dTheta_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( theta, stress );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dTheta_dStress(
                                           stress )
                                           .norm(),
                                         dTheta_dStress_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dJ2_dStress()
{

  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };

  // lambda function to calculate J2
  const auto J2 = []( const Vector6d& stress ) {
    Eigen::MatrixXd j2( 1, 1 );
    j2 << ( Marmot::ContinuumMechanics::VoigtNotation::Invariants::J2( stress ) );
    return j2;
  };

  const auto dJ2_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( J2, stress );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dJ2_dStress( stress )
                                           .norm(),
                                         dJ2_dStress_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dJ3_dStress()
{

  const Vector6d stress = { 1, 2, 3, 4, 5, 6 };

  // lambda function to calculate J3
  const auto J3 = []( const Vector6d& stress ) {
    Eigen::MatrixXd j3( 1, 1 );
    j3 << ( Marmot::ContinuumMechanics::VoigtNotation::Invariants::J3( stress ) );
    return j3;
  };

  const auto dJ3_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( J3, stress );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dJ3_dStress( stress )
                                           .norm(),
                                         dJ3_dStress_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dJ2Strain_dStrain()
{

  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };

  // lambda function to calculate J2
  const auto J2 = []( const Vector6d& strain ) {
    Eigen::MatrixXd j2( 1, 1 );
    j2 << ( Marmot::ContinuumMechanics::VoigtNotation::Invariants::J2Strain( strain ) );
    return j2;
  };

  const auto dJ2_dStrain_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( J2, strain );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dJ2Strain_dStrain(
                                           strain )
                                           .norm(),
                                         dJ2_dStrain_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dJ3Strain_dStrain()
{
  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };

  // lambda function to calculate J3
  const auto J3 = []( const Vector6d& strain ) {
    Eigen::MatrixXd j3( 1, 1 );
    j3 << ( Marmot::ContinuumMechanics::VoigtNotation::Invariants::J3Strain( strain ) );
    return j3;
  };

  const auto dJ3_dStrain_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( J3, strain );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::dJ3Strain_dStrain(
                                           strain )
                                           .norm(),
                                         dJ3_dStrain_FD.norm(),
                                         1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void test_dSortedPrincipalStrains_dStrain()
{

  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };

  // lambda function to calculate sorted principal strains
  const auto sortedPrincipalStrains = []( const Vector6d& strain ) {
    Eigen::MatrixXd
      sortedPrincipalStrains = Marmot::ContinuumMechanics::VoigtNotation::Invariants::sortedPrincipalStrains( strain );
    return sortedPrincipalStrains;
  };

  const auto dSortedPrincipalStrains_dStrain_FD = Marmot::NumericalAlgorithms::Differentiation::
    forwardDifference( sortedPrincipalStrains, strain );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::
                                           dSortedStrainPrincipal_dStrain( strain )
                                             .norm(),
                                         dSortedPrincipalStrains_dStrain_FD.norm(),
                                         1e-5 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessToVoigtRoundTripEigenTensor()
{
  Marmot::Matrix6d voigtStiffness;
  // clang-format off
  voigtStiffness <<
    1200,   400,    400,    50,     60,     70,
    400,    1200,   400,    80,     90,     100,
    400,    400,    1200,   110,    120,    130,
    50,     80,     110,    400,    140,    150,
    60,     90,     120,    140,    400,    160,
    70,     100,    130,    150,    160,    400;
  // clang-format on

  const auto
    stiffnessTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness< Marmot::EigenTensors::Tensor3333d >(
      voigtStiffness );
  const auto voigtStiffnessRoundTrip = Marmot::ContinuumMechanics::VoigtNotation::stiffnessToVoigt( stiffnessTensor );

  throwExceptionOnFailure( checkIfEqual< double >( voigtStiffnessRoundTrip, voigtStiffness, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessToVoigtRoundTripFastorTensor()
{
  Marmot::Matrix6d voigtStiffnessEigen;
  // clang-format off
  voigtStiffnessEigen <<
    1200,   400,    400,    50,     60,     70,
    400,    1200,   400,    80,     90,     100,
    400,    400,    1200,   110,    120,    130,
    50,     80,     110,    400,    140,    150,
    60,     90,     120,    140,    400,    160,
    70,     100,    130,    150,    160,    400;
  // clang-format on

  const Fastor::Tensor< double, 6, 6 > voigtStiffness( voigtStiffnessEigen.data(), Fastor::ColumnMajor );

  const auto stiffnessTensor         = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( voigtStiffness );
  const auto voigtStiffnessRoundTrip = Marmot::ContinuumMechanics::VoigtNotation::stiffnessToVoigt( stiffnessTensor );

  throwExceptionOnFailure( checkIfEqual< double >( voigtStiffnessRoundTrip, voigtStiffnessEigen, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStiffnessToVoigtAveragesMinorSymmetricPermutations()
{
  // Build a fourth-order tensor directly (NOT via voigtToStiffness, which always produces an
  // already minor-symmetric tensor) with deliberately broken minor symmetry, so that
  // stiffnessToVoigt's averaging of the four minor-symmetric-equivalent entries is actually
  // exercised.
  Marmot::EigenTensors::Tensor3333d C;
  C.setZero();

  // Voigt entry (a,b) = (3,3), i.e. (i,j) = (0,1), (k,l) = (0,1): all four permutations distinct.
  C( 0, 1, 0, 1 ) = 10;
  C( 1, 0, 0, 1 ) = 20;
  C( 1, 0, 1, 0 ) = 30;
  C( 0, 1, 1, 0 ) = 40;
  // expected voigtStiffness(3,3) = (10 + 20 + 30 + 40) / 4 = 25

  // Voigt entry (a,b) = (0,3), i.e. (i,j) = (0,0), (k,l) = (0,1): the (i,j) swap term is
  // degenerate (i == j), so only the (k,l) swap contributes a second distinct value.
  C( 0, 0, 0, 1 ) = 100;
  C( 0, 0, 1, 0 ) = 200;
  // expected voigtStiffness(0,3) = (100 + 100 + 200 + 200) / 4 = 150

  const auto voigtStiffness = Marmot::ContinuumMechanics::VoigtNotation::stiffnessToVoigt( C );

  throwExceptionOnFailure( checkIfEqual( voigtStiffness( 3, 3 ), 25.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed at (3,3)" );
  throwExceptionOnFailure( checkIfEqual( voigtStiffness( 0, 3 ), 150.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed at (0,3)" );
}

void testVoigtToStiffnessMatchesAnalyticIsotropicTensor()
{
  // Independent ground truth: the classical isotropic stiffness tensor
  // C_ijkl = lambda * delta_ij * delta_kl + mu * ( delta_ik * delta_jl + delta_il * delta_jk ),
  // built here from a plain second-order identity tensor via Fastor tensor algebra (not reusing
  // any Marmot tensor-building helper), and compared against
  // voigtToStiffness( Isotropic::stiffnessTensor( E, nu ) ).
  const double E  = 1000.;
  const double nu = 0.25;

  const double lambda = nu * E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
  const double mu     = E / ( 2 * ( 1 + nu ) );

  Fastor::Tensor< double, 3, 3 > delta;
  delta.eye2();

  using namespace FastorIndices;
  FastorStandardTensors::Tensor3333d expectedStiffnessTensor = lambda * Fastor::outer( delta, delta ) +
                                                               mu *
                                                                 ( Fastor::einsum< ik, jl, to_ijkl >( delta, delta ) +
                                                                   Fastor::einsum< il, jk, to_ijkl >( delta, delta ) );

  const Marmot::Matrix6d voigtStiffness = Marmot::ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

  // voigtToStiffness defaults to a Fastor result, so this also exercises the Eigen-in/Fastor-out
  // (cross-library) conversion path.
  const auto stiffnessTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( voigtStiffness );

  throwExceptionOnFailure( checkIfEqual( stiffnessTensor, expectedStiffnessTensor, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testVoigtToStiffnessIndexMapping()
{
  // Spot-check that voigtToStiffness places each Voigt entry at the tensor indices implied by the
  // standard Voigt ordering (0->xx, 1->yy, 2->zz, 3->xy, 4->xz, 5->yz), with a focus on the
  // shear-index mapping (12/13/23 <-> 3/4/5), which is the most error-prone part of the mapping.
  Marmot::Matrix6d voigtStiffness;
  // clang-format off
  voigtStiffness <<
    1200,   400,    400,    50,     60,     70,
    400,    1200,   400,    80,     90,     100,
    400,    400,    1200,   110,    120,    130,
    50,     80,     110,    400,    140,    150,
    60,     90,     120,    140,    400,    160,
    70,     100,    130,    150,    160,    400;
  // clang-format on

  const auto
    stiffnessTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness< Marmot::EigenTensors::Tensor3333d >(
      voigtStiffness );

  throwExceptionOnFailure( checkIfEqual( stiffnessTensor( 0, 1, 0, 1 ), voigtStiffness( 3, 3 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for (xy,xy) -> (3,3)" );
  throwExceptionOnFailure( checkIfEqual( stiffnessTensor( 2, 0, 1, 2 ), voigtStiffness( 4, 5 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for (xz,yz) -> (4,5)" );
  throwExceptionOnFailure( checkIfEqual( stiffnessTensor( 1, 2, 2, 0 ), voigtStiffness( 5, 4 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for (yz,xz) -> (5,4)" );
  throwExceptionOnFailure( checkIfEqual( stiffnessTensor( 0, 0, 1, 2 ), voigtStiffness( 0, 5 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for (xx,yz) -> (0,5)" );
  throwExceptionOnFailure( checkIfEqual( stiffnessTensor( 1, 1, 2, 0 ), voigtStiffness( 1, 4 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for (yy,xz) -> (1,4)" );
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple, previously-untested conversions and invariants
// ─────────────────────────────────────────────────────────────────────────────

void testVoigtToAxisymmetricVoigt()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  Vector6d voigt = { 1, 2, 3, 4, 5, 6 };
  throwExceptionOnFailure( checkIfEqual< double >( voigtToAxisymmetricVoigt( voigt ), Vector4d( 1, 2, 3, 4 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testAxisymmetricVoigtToVoigt()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  Vector4d voigtAxisym = { 1, 2, 3, 4 };
  throwExceptionOnFailure( checkIfEqual< double >( axisymmetricVoigtToVoigt( voigtAxisym ),
                                                   Vector6d( 1, 2, 3, 4, 0, 0 ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testNormStress()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // A pure-shear Voigt stress (3,4) maps to a 3x3 stress matrix with off-diagonal 3 and 4 (each
  // appearing twice, symmetric), so ||sigma||_F = sqrt(2*3^2 + 2*4^2).
  Vector6d stress = { 0, 0, 0, 3, 4, 0 };
  throwExceptionOnFailure( checkIfEqual( Invariants::normStress( stress ), std::sqrt( 2. * 9. + 2. * 16. ), 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testStrainVolumetricNegative()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // A diagonal (already-principal) strain state with two negative and one positive principal
  // strain: the negative part is macauly(-(-0.01)) + macauly(-(-0.02)) + macauly(-(0.03)) = 0.01+0.02+0.
  Vector6d strain = { -0.01, -0.02, 0.03, 0, 0, 0 };
  throwExceptionOnFailure( checkIfEqual( Invariants::StrainVolumetricNegative( strain ), 0.03, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testI1Strain()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  Vector6d strain = { 1, 2, 3, 4, 5, 6 };
  throwExceptionOnFailure( checkIfEqual( Invariants::I1Strain( strain ), 6.0, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

// ─────────────────────────────────────────────────────────────────────────────
// Transformations namespace: entirely untested previously. Verified via physical/algebraic
// identities (Cauchy's theorem for the projection matrices, round-trip local<->global, rotational
// invariance of an isotropic stiffness) rather than re-deriving the transformation formulas.
// ─────────────────────────────────────────────────────────────────────────────

namespace {
  // A fixed, non-trivial (non-axis-aligned) orthonormal coordinate system, used throughout the
  // Transformations tests below.
  Eigen::Matrix3d makeTestRotation()
  {
    Eigen::Matrix3d N;
    // clang-format off
    N << 0.7071067811865476,  0.7071067811865476, 0.0,
        -0.5,                 0.5,                 0.7071067811865476,
         0.5,                -0.5,                 0.7071067811865476;
    // clang-format on
    return N;
  }
} // namespace

void testTransformationMatrixStressVoigtIsIdentityForIdentitySystem()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Matrix6d T = Transformations::transformationMatrixStressVoigt( Eigen::Matrix3d::Identity() );
  throwExceptionOnFailure( checkIfEqual< double >( T, Matrix6d::Identity() ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testTransformationMatrixStressVoigtMatchesRotateVoigtStress()
{
  // transformationMatrixStressVoigt(X) internally uses directionCosines(X) = X^T as the
  // rows-are-new-basis-vectors matrix, i.e. it is equivalent to rotateVoigtStress(X^T, ...).
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Eigen::Matrix3d N      = makeTestRotation();
  const Vector6d        stress = { 10, -20, 30, 4, -5, 6 };

  const Matrix6d T           = Transformations::transformationMatrixStressVoigt( N );
  const Vector6d viaMatrix   = T * stress;
  const Vector6d viaRotation = Transformations::rotateVoigtStress( N.transpose(), stress );

  throwExceptionOnFailure( checkIfEqual< double >( viaMatrix, viaRotation, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " transformationMatrixStressVoigt() does not match rotateVoigtStress()" );
}

void testTransformationMatrixStrainVoigtMatchesDirectStrainRotation()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Eigen::Matrix3d N      = makeTestRotation();
  const Vector6d        strain = { 0.001, -0.002, 0.003, 0.0004, -0.0005, 0.0006 };

  const Matrix6d T         = Transformations::transformationMatrixStrainVoigt( N );
  const Vector6d viaMatrix = T * strain;

  // See testTransformationMatrixStressVoigtMatchesRotateVoigtStress(): the matrix actually applied
  // is directionCosines(N) = N^T.
  const Matrix3d strainMat  = voigtToStrain( strain );
  const Matrix3d rotatedMat = N.transpose() * strainMat * N;
  const Vector6d viaDirect  = strainToVoigt( rotatedMat );

  throwExceptionOnFailure( checkIfEqual< double >( viaMatrix, viaDirect, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " transformationMatrixStrainVoigt() does not match N^T*strain*N" );
}

void testProjectVoigtStressToPlaneMatchesCauchyTraction()
{
  // Cauchy's stress theorem: t = sigma * n. projectVoigtStressToPlane(n) * stressVoigt must give
  // the same traction vector as the direct matrix-vector product.
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector6d stress    = { 10, -20, 30, 4, -5, 6 };
  const Matrix3d stressMat = voigtToStress( stress );
  const Vector3d n         = Vector3d( 1., 2., -2. ).normalized();

  const Vector3d tractionDirect        = stressMat * n;
  const Vector3d tractionViaProjection = Transformations::projectVoigtStressToPlane( n ) * stress;

  throwExceptionOnFailure( checkIfEqual< double >( tractionViaProjection, tractionDirect, 1e-10 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " projectVoigtStressToPlane() does not reproduce t = sigma*n" );
}

void testProjectVoigtStrainToPlaneMatchesDirectProduct()
{
  // There is no equally simple physical identity for the strain projection (the factor-of-2
  // engineering-shear convention breaks the clean t=sigma*n analogy), so this only checks the
  // documented relation to projectVoigtStressToPlane(): the off-diagonal (top-right) block is
  // halved, everything else identical.
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector3d  n          = Vector3d( 1., 2., -2. ).normalized();
  const Matrix36d stressProj = Transformations::projectVoigtStressToPlane( n );
  const Matrix36d strainProj = Transformations::projectVoigtStrainToPlane( n );

  Matrix36d expected = stressProj;
  expected.topRightCorner( 3, 3 ) *= 0.5;

  throwExceptionOnFailure( checkIfEqual< double >( strainProj, expected, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testRotateVoigtStressRoundTrip()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Eigen::Matrix3d N      = makeTestRotation();
  const Vector6d        stress = { 10, -20, 30, 4, -5, 6 };

  const Vector6d rotated     = Transformations::rotateVoigtStress( N, stress );
  const Vector6d rotatedBack = Transformations::rotateVoigtStress( N.transpose(), rotated );

  throwExceptionOnFailure( checkIfEqual< double >( rotatedBack, stress, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " rotating by N then N^-1 did not recover the original stress" );
}

void testTransformStressStrainLocalGlobalRoundTrip()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Eigen::Matrix3d N      = makeTestRotation();
  const Vector6d        stress = { 10, -20, 30, 4, -5, 6 };
  const Vector6d        strain = { 0.001, -0.002, 0.003, 0.0004, -0.0005, 0.0006 };

  const Vector6d stressLocal     = Transformations::transformStressToLocalSystem( stress, N );
  const Vector6d stressRoundTrip = Transformations::transformStressToGlobalSystem( stressLocal, N );
  throwExceptionOnFailure( checkIfEqual< double >( stressRoundTrip, stress, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " stress local/global round-trip failed" );

  const Vector6d strainLocal     = Transformations::transformStrainToLocalSystem( strain, N );
  const Vector6d strainRoundTrip = Transformations::transformStrainToGlobalSystem( strainLocal, N );
  throwExceptionOnFailure( checkIfEqual< double >( strainRoundTrip, strain, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " strain local/global round-trip failed" );
}

void testTransformStiffnessToGlobalSystemPreservesIsotropicStiffness()
{
  // An isotropic stiffness tensor must be invariant under any rotation.
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Matrix6d        C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 20000., 0.25 );
  const Eigen::Matrix3d N = makeTestRotation();

  const Matrix6d rotated = Transformations::transformStiffnessToGlobalSystem( C, N );

  throwExceptionOnFailure( checkIfEqual< double >( rotated, C, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " an isotropic stiffness must be invariant under rotation" );
}

// ─────────────────────────────────────────────────────────────────────────────
// Invariants::principalValuesAndDerivatives(): a "fast implementation of the classical algorithm"
// for the eigenvalues (and their derivatives) of a symmetric 3x3 matrix in Voigt notation, with no
// callers anywhere in the codebase. Verified against Eigen::SelfAdjointEigenSolver (already
// trusted, used by principalStresses() above) for the values, and against numerical
// differentiation of the function's own value output for the derivatives.
// ─────────────────────────────────────────────────────────────────────────────

void testPrincipalValuesAndDerivativesGeneralCase()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector6d stress = { 10., -20., 30., 4., -5., 6. };

  const auto [e, dE_dS] = Invariants::principalValuesAndDerivatives( stress );

  // Values: compare against the sorted eigenvalues from the already-trusted SelfAdjointEigenSolver
  // path (order-independent, since e's own ordering is a formula artifact, not sorted).
  Vector3d eSorted = e;
  std::sort( eSorted.data(), eSorted.data() + 3 );
  Vector3d trueSorted = Invariants::principalStresses( stress ); // SelfAdjointEigenSolver eigenvalues, ascending
  throwExceptionOnFailure( checkIfEqual< double >( eSorted, trueSorted, 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " principal values do not match Eigen::SelfAdjointEigenSolver" );

  // Derivatives: central-difference the function's own value output w.r.t. each Voigt component.
  const double h = 1e-6;
  for ( int j = 0; j < 6; j++ ) {
    Vector6d sp = stress, sm = stress;
    sp( j ) += h;
    sm( j ) -= h;
    const auto [ep, dEp_dS] = Invariants::principalValuesAndDerivatives( sp );
    const auto [em, dEm_dS] = Invariants::principalValuesAndDerivatives( sm );
    (void)dEp_dS;
    (void)dEm_dS;
    const Vector3d dE_dSj_num = ( ep - em ) / ( 2. * h );

    throwExceptionOnFailure( checkIfEqual< double >( dE_dS.col( j ), dE_dSj_num, 1e-4 ),
                             MakeString() << __PRETTY_FUNCTION__ << " column " << j
                                          << " does not match the numerical derivative" );
  }
}

void testPrincipalValuesAndDerivativesDiagonalStress()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // p1 = S(3)^2+S(4)^2+S(5)^2 == 0: hits the "matrix is already diagonal" early return.
  const Vector6d stress = { 1., 2., 3., 0., 0., 0. };

  const auto [e, dE_dS] = Invariants::principalValuesAndDerivatives( stress );

  throwExceptionOnFailure( checkIfEqual< double >( e, Vector3d( 1., 2., 3. ), 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " must return the diagonal entries directly, in order" );

  Eigen::Matrix< double, 3, 6 > dE_dS_expected;
  // clang-format off
  dE_dS_expected << 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0;
  // clang-format on
  throwExceptionOnFailure( checkIfEqual< double >( dE_dS, dE_dS_expected, 1e-14 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed for the diagonal-input derivative" );
}

void testPrincipalValuesAndDerivativesTriaxialRGreaterEqualOne()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // diag(1,1,3) rotated 45 degrees about the y-axis: same eigenvalues {1,1,3} as the diagonal
  // input, but with a nonzero off-diagonal S13 term (so p1 > 0, avoiding the diagonal shortcut
  // above) chosen so that, by construction, r == +1 exactly -- the "r >= 1" boundary branch.
  const Vector6d stress = { 2., 1., 2., 0., 1., 0. };

  const auto [e, dE_dS] = Invariants::principalValuesAndDerivatives( stress );
  (void)dE_dS; // dPhi_dR is deliberately clamped to 0 at this boundary, not a smooth limit -- only
               // the values are checked here, not the derivative.

  Vector3d eSorted = e;
  std::sort( eSorted.data(), eSorted.data() + 3 );
  throwExceptionOnFailure( checkIfEqual< double >( eSorted, Vector3d( 1., 1., 3. ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__ << " principal values are wrong at the r>=1 boundary" );
}

void testPrincipalValuesAndDerivativesTriaxialNearRMinusOne()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // diag(3,3,1) rotated 45 degrees about the y-axis: eigenvalues {3,3,1}, constructed so that r is
  // mathematically exactly -1. In practice, floating-point rounding through acos()/sqrt() in the
  // implementation lands the computed r only ~1e-7 short of -1 for every construction tried, not
  // within the "r <= -1" branch's guard band -- unlike the analogous r>=1 case, which reliably
  // rounds past its boundary. So this test only exercises the smooth interior acos(r)/3 formula in
  // its near-singular regime, not the dedicated "r <= -1" branch itself; it checks that the
  // principal values are still correct there, not that the boundary branch ran.
  const Vector6d stress = { 2., 3., 2., 0., -1., 0. };

  const auto [e, dE_dS] = Invariants::principalValuesAndDerivatives( stress );
  (void)dE_dS;

  Vector3d eSorted = e;
  std::sort( eSorted.data(), eSorted.data() + 3 );
  throwExceptionOnFailure( checkIfEqual< double >( eSorted, Vector3d( 1., 3., 3. ), 1e-8 ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " principal values are wrong near the r<=-1 boundary" );
}

// ─────────────────────────────────────────────────────────────────────────────
// The Lode angle derivatives are singular at the triaxial boundaries (theta == 0 or theta == Pi/3),
// where sin(3*theta) == 0. The implementations regularize this case; all that is required there is a
// finite result, so that no inf/NaN propagates into a return mapping or tangent. The exactly-triaxial
// states below (principal values {1,1,3}) make haighWestergaard()/haighWestergaardFromStrain() assign
// theta the exact literal 0.
// ─────────────────────────────────────────────────────────────────────────────

void testDThetaDStressIsFiniteAtLodeAngleBoundary()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector6d stressThetaZero = { 2., 1., 2., 0., 1., 0. };

  throwExceptionOnFailure( ContinuumMechanics::HaighWestergaard::haighWestergaard( stressThetaZero ).theta == 0.0,
                           MakeString() << __PRETTY_FUNCTION__ << " test stress does not have theta == 0" );

  throwExceptionOnFailure( std::isfinite( Derivatives::dTheta_dJ2( stressThetaZero ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " dTheta_dJ2 not finite for theta == 0" );
  throwExceptionOnFailure( std::isfinite( Derivatives::dTheta_dJ3( stressThetaZero ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " dTheta_dJ3 not finite for theta == 0" );
  throwExceptionOnFailure( Derivatives::dTheta_dStress( stressThetaZero ).allFinite(),
                           MakeString() << __PRETTY_FUNCTION__ << " dTheta_dStress not finite for theta == 0" );
}

void testDThetaStrainDStrainIsFiniteAtLodeAngleBoundary()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // Same underlying tensor as the stress case above, but as a Voigt strain vector (engineering shear).
  const Vector6d strainThetaZero = { 2., 1., 2., 0., 2., 0. };

  throwExceptionOnFailure( ContinuumMechanics::HaighWestergaard::haighWestergaardFromStrain( strainThetaZero ).theta ==
                             0.0,
                           MakeString() << __PRETTY_FUNCTION__ << " test strain does not have theta == 0" );

  throwExceptionOnFailure( std::isfinite( Derivatives::dThetaStrain_dJ2Strain( strainThetaZero ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " dThetaStrain_dJ2Strain not finite for theta == 0" );
  throwExceptionOnFailure( std::isfinite( Derivatives::dThetaStrain_dJ3Strain( strainThetaZero ) ),
                           MakeString() << __PRETTY_FUNCTION__ << " dThetaStrain_dJ3Strain not finite for theta == 0" );
  throwExceptionOnFailure( Derivatives::dThetaStrain_dStrain( strainThetaZero ).allFinite(),
                           MakeString() << __PRETTY_FUNCTION__ << " dThetaStrain_dStrain not finite for theta == 0" );
}

// ─────────────────────────────────────────────────────────────────────────────
// The remaining Derivatives-namespace functions below have no callers anywhere in the codebase
// (they appear to be infrastructure ported for damage-plasticity models not present in this
// open-source build, per their "PhD Thesis David Unteregger" / "Code David" references).
// ─────────────────────────────────────────────────────────────────────────────

void testDThetaStrainDStrainMatchesNumericalDifferentiation()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation::Derivatives;
  const Vector6d strain = { 1, 2, 3, 4, 5, 6 };

  const auto theta = []( const Vector6d& strain ) {
    Eigen::MatrixXd theta( 1, 1 );
    theta << ContinuumMechanics::HaighWestergaard::haighWestergaardFromStrain( strain ).theta;
    return theta;
  };

  const auto dTheta_dStrain_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( theta, strain );

  throwExceptionOnFailure( checkIfEqual< double >( dThetaStrain_dStrain( strain ).transpose(),
                                                   dTheta_dStrain_FD,
                                                   1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testDStressPrincipalsDStressMatchesIndependentCentralDifference()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector6d stress = { 10., -20., 30., 4., -5., 6. };

  const Matrix36 J = Derivatives::dStressPrincipals_dStress( stress );

  // Independently re-derive the same central-difference Jacobian with a different (fixed) step
  // size, rather than trusting the function's own choice of step.
  const double h = 1e-5;
  Matrix36     J_independent;
  for ( int i = 0; i < 6; i++ ) {
    Vector6d left = stress, right = stress;
    left( i ) -= h;
    right( i ) += h;
    J_independent.col( i ) = ( Invariants::principalStresses( right ) - Invariants::principalStresses( left ) ) /
                             ( 2. * h );
  }

  throwExceptionOnFailure( checkIfEqual< double >( J, J_independent, 1e-4 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

void testDStrainVolumetricNegativeDStrainPrincipal()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // A diagonal strain with one negative and two positive principal strains.
  const Vector6d strain = { -0.01, 0.02, 0.03, 0, 0, 0 };

  const Vector3d dEvdEpPrinc = Derivatives::dStrainVolumetricNegative_dStrainPrincipal( strain );

  // Cross-check against the actual sorted principal strains used internally, rather than assuming
  // their order matches the input.
  const Vector3d sorted = Invariants::sortedPrincipalStrains( strain );
  for ( int i = 0; i < 3; i++ ) {
    const double expected = sorted( i ) < 0.0 ? -1.0 : 0.0;
    throwExceptionOnFailure( checkIfEqual( dEvdEpPrinc( i ), expected, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << " failed at index " << i );
  }
}

void testDEpDEAndDDeltaEpvDE()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // CelInv * Cep == Identity: purely elastic increment, no plastic strain accrues.
  throwExceptionOnFailure( checkIfEqual< double >( Derivatives::dEp_dE( Matrix6d::Identity(), Matrix6d::Identity() ),
                                                   Matrix6d::Zero() ),
                           MakeString() << __PRETTY_FUNCTION__ << " dEp_dE failed for CelInv*Cep == Identity" );

  // CelInv * Cep == Zero: fully plastic increment (Cep == 0), so dEp/dE == Identity.
  throwExceptionOnFailure( checkIfEqual< double >( Derivatives::dEp_dE( Matrix6d::Identity(), Matrix6d::Zero() ),
                                                   Matrix6d::Identity() ),
                           MakeString() << __PRETTY_FUNCTION__ << " dEp_dE failed for Cep == 0" );

  // dDeltaEpv_dE must equal I^T * dEp_dE by definition.
  const Matrix6d CelInv = ( Matrix6d() << 1,
                            0.1,
                            0,
                            0,
                            0,
                            0,
                            0.1,
                            1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1 )
                            .finished();
  const Matrix6d Cep = 0.5 * Matrix6d::Identity();
  throwExceptionOnFailure( checkIfEqual< double >( Derivatives::dDeltaEpv_dE( CelInv, Cep ),
                                                   ( I.transpose() * Derivatives::dEp_dE( CelInv, Cep ) ).eval(),
                                                   1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " dDeltaEpv_dE failed" );
}

void testDSortedStrainPrincipalDStrainHydrostaticBranch()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  // A purely hydrostatic strain: rho == 0, hitting the "near-origin" fallback branch (the 1e16
  // sentinel documented as ported from "Code David").
  const Vector6d strain = { 0.01, 0.01, 0.01, 0, 0, 0 };

  throwExceptionOnFailure( ContinuumMechanics::HaighWestergaard::haighWestergaardFromStrain( strain ).rho == 0.0,
                           MakeString() << __PRETTY_FUNCTION__ << " test strain is not exactly hydrostatic" );

  const Matrix36 dEpPrincdEp = Derivatives::dSortedStrainPrincipal_dStrain( strain );

  // dEpPrinc_dEprho contributes sqrt(2/3)*cos(...) (bounded) times the 1e16 sentinel row, so every
  // entry of the result must be very large in magnitude.
  throwExceptionOnFailure( dEpPrincdEp.cwiseAbs().minCoeff() > 1e10,
                           MakeString() << __PRETTY_FUNCTION__
                                        << " expected very large entries at the hydrostatic sentinel branch" );
}

void testDDeltaEpvnegDEIsConsistentWithItsComponents()
{
  using namespace Marmot::ContinuumMechanics::VoigtNotation;
  const Vector6d strain = { -0.01, 0.02, 0.03, 0.004, -0.005, 0.006 };
  const Matrix6d CelInv = Matrix6d::Identity();
  const Matrix6d Cep    = 0.3 * Matrix6d::Identity();

  const RowVector6d result = Derivatives::dDeltaEpvneg_dE( strain, CelInv, Cep );

  const RowVector6d expected = Derivatives::dStrainVolumetricNegative_dStrainPrincipal( strain ).transpose() *
                               Derivatives::dSortedStrainPrincipal_dStrain( strain ) *
                               Derivatives::dEp_dE( CelInv, Cep );

  throwExceptionOnFailure( checkIfEqual< double >( result, expected, 1e-12 ),
                           MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ testStrainToVoigt,
                                                       testVoigtToPlaneVoigt,
                                                       testPrincipalStrains,
                                                       testPrincipalStresses,
                                                       testSortedPrincipalStrains,
                                                       testPrincipalStressDirections,
                                                       testVonMisesEquivalentStress,
                                                       testVonMisesEquivalentStrain,
                                                       testI1,
                                                       testI2,
                                                       testI2Strain,
                                                       testI3,
                                                       testI3Strain,
                                                       testJ2,
                                                       testJ3,
                                                       testJ3Strain,
                                                       test_dStressMean_dStress,
                                                       test_dRho_dStress,
                                                       test_dRhoStrain_dStrain,
                                                       test_dTheta_dStress,
                                                       test_dJ2_dStress,
                                                       test_dJ3_dStress,
                                                       test_dJ2Strain_dStrain,
                                                       test_dJ3Strain_dStrain,
                                                       test_dSortedPrincipalStrains_dStrain,
                                                       testStiffnessToVoigtRoundTripEigenTensor,
                                                       testStiffnessToVoigtRoundTripFastorTensor,
                                                       testStiffnessToVoigtAveragesMinorSymmetricPermutations,
                                                       testVoigtToStiffnessMatchesAnalyticIsotropicTensor,
                                                       testVoigtToStiffnessIndexMapping,
                                                       testVoigtToAxisymmetricVoigt,
                                                       testAxisymmetricVoigtToVoigt,
                                                       testNormStress,
                                                       testStrainVolumetricNegative,
                                                       testI1Strain,
                                                       testTransformationMatrixStressVoigtIsIdentityForIdentitySystem,
                                                       testTransformationMatrixStressVoigtMatchesRotateVoigtStress,
                                                       testTransformationMatrixStrainVoigtMatchesDirectStrainRotation,
                                                       testProjectVoigtStressToPlaneMatchesCauchyTraction,
                                                       testProjectVoigtStrainToPlaneMatchesDirectProduct,
                                                       testRotateVoigtStressRoundTrip,
                                                       testTransformStressStrainLocalGlobalRoundTrip,
                                                       testTransformStiffnessToGlobalSystemPreservesIsotropicStiffness,
                                                       testPrincipalValuesAndDerivativesGeneralCase,
                                                       testPrincipalValuesAndDerivativesDiagonalStress,
                                                       testPrincipalValuesAndDerivativesTriaxialRGreaterEqualOne,
                                                       testPrincipalValuesAndDerivativesTriaxialNearRMinusOne,
                                                       testDThetaDStressIsFiniteAtLodeAngleBoundary,
                                                       testDThetaStrainDStrainIsFiniteAtLodeAngleBoundary,
                                                       testDThetaStrainDStrainMatchesNumericalDifferentiation,
                                                       testDStressPrincipalsDStressMatchesIndependentCentralDifference,
                                                       testDStrainVolumetricNegativeDStrainPrincipal,
                                                       testDEpDEAndDDeltaEpvDE,
                                                       testDSortedStrainPrincipalDStrainHydrostaticBranch,
                                                       testDDeltaEpvnegDEIsConsistentWithItsComponents };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
