#include "Marmot/HaighWestergaard.h"
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
                                         dJ2_dStress_FD.norm() ),
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
                                         dJ3_dStress_FD.norm() ),
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
                                         dJ2_dStrain_FD.norm() ),
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
                                         dJ3_dStrain_FD.norm() ),
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

int main()
{
  testStrainToVoigt();
  testVoigtToPlaneVoigt();
  testPrincipalStrains();
  testPrincipalStresses();
  testSortedPrincipalStrains();
  testPrincipalStressDirections();
  testVonMisesEquivalentStress();
  testVonMisesEquivalentStrain();
  testI1();
  testI2();
  testI2Strain();
  testI3();
  testI3Strain();
  testJ2();
  testJ3();
  testJ3Strain();

  // Derivatives test
  test_dStressMean_dStress();
  test_dRho_dStress();
  test_dTheta_dStress();
  test_dJ2_dStress();
  test_dJ3_dStress();
  test_dJ2Strain_dStrain();
  test_dJ3Strain_dStrain();

  test_dSortedPrincipalStrains_dStrain();

  return 0;
}
