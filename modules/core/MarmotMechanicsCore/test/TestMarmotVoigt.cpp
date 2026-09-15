#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
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
  const auto     hw     = ContinuumMechanics::HaighWestergaard::haighWestergaard( stress );

  const auto theta = []( const Vector6d& stress ) {
    Eigen::MatrixXd theta( 1, 1 );
    theta << ContinuumMechanics::HaighWestergaard::haighWestergaard( stress ).theta;
    return theta;
  };

  const auto dTheta_dStress_FD = Marmot::NumericalAlgorithms::Differentiation::forwardDifference( theta, stress );

  throwExceptionOnFailure( checkIfEqual( Marmot::ContinuumMechanics::VoigtNotation::Derivatives::
                                           dTheta_dStress( hw.theta, stress )
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
                                                       testVoigtToStiffnessIndexMapping };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
