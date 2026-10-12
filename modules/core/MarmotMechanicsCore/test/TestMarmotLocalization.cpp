#include "Marmot/MarmotLocalization.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/VonMises.h"

using namespace Marmot::Testing;

void testMarmotLocalization()
{
  using namespace Marmot::ContinuumMechanics::LocalizationAnalysis;

  // material properties
  Eigen::Vector< double, 6 > materialProperties;
  materialProperties << 210000., 0.3, 200., 0., 0., 0.;

  auto material = Marmot::Materials::VonMisesModel( materialProperties.data(), materialProperties.size(), 0 );

  // set state vars
  if ( material.getNumberOfRequiredStateVars() > 1 ) {
    throw std::runtime_error( "Number of required state vars for Mises model changed!" );
  }
  double kappa = 0;

  Marmot::Vector6d stress;
  Marmot::Matrix6d dStress_dStrain;
  Marmot::Vector6d dStrain;
  const double     dT = 1;

  // initialize stress on shear meridian close to yield surface
  stress << 115.4, -115.4, 0., 0., 0., 0.;

  // set strain increment to reach shear meridian
  dStrain << 1e-6, -1e-6, 0., 0., 0., 0.;

  // set material state
  MarmotMaterialHypoElastic::state3D state;
  state.stress               = stress;
  state.elasticEnergyDensity = 0.;
  state.dissipation          = 0.;
  state.stateVars            = &kappa;

  // set time info
  MarmotMaterialHypoElastic::timeInfo timeInfo;
  timeInfo.time = 1.0;
  timeInfo.dT   = dT;
  // evaluate material
  material.computeStress( state, dStress_dStrain, dStrain, timeInfo );

  // normal vector
  Marmot::Vector3d n;
  n << 1, -1, 0;
  n /= n.norm();

  throwExceptionOnFailure( checkIfEqual( computeAcousticTensor( dStress_dStrain, n ).determinant(), 0. ),
                           "for n = [1/√2, 1/√2, 0], determinant of acoustic tensor should be zero on shear meridian" );

  // normal vector
  n << 1, 0, 0;
  n /= n.norm();

  throwExceptionOnFailure( checkIfEqual( minimumDeterminantAcousticTensor( dStress_dStrain ) /
                                           computeAcousticTensor( dStress_dStrain, n ).determinant(),
                                         0.0301537,
                                         1e-6 ),
                           "test minimumDeterminantAcousticTensor" );
}

void testLocalizationChecker()
{
  using namespace Marmot::ContinuumMechanics::LocalizationAnalysis;

  // A well-conditioned (positive-definite, non-singular) acoustic tensor: no localization.
  Marmot::Matrix3d wellConditioned = Marmot::Matrix3d::Identity() * 1000.;
  throwExceptionOnFailure( localizationChecker( wellConditioned ) == false,
                           "localizationChecker() must return false for a well-conditioned acoustic tensor in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // A singular (zero-determinant) acoustic tensor: localization detected.
  Marmot::Matrix3d singular = Marmot::Matrix3d::Zero();
  singular( 0, 0 )          = 1000.;
  singular( 1, 1 )          = 1000.;
  // row/col 2 left as zero -> determinant exactly 0
  throwExceptionOnFailure( localizationChecker( singular ) == true,
                           "localizationChecker() must return true for a singular acoustic tensor in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testMarmotLocalization, testLocalizationChecker };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
