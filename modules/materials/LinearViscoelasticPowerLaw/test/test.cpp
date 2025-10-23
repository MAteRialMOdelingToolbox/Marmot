#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

void testLinearViscoelasticPowerLaw()
{
  // material properties
  const int LVPLCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( "LINEARVISCOELASTICPOWERLAW" );

  // material properties
  Eigen::Vector< double, 7 > materialProperties;
  // elastic parameters
  double E  = 2e5;
  double nu = 0.2;
  // viscoelastic parameters
  double m          = 0.5;
  double n          = 0.1;
  int    nKelvin    = 10;
  double minTau     = 0.0001;
  double timeToDays = 1.;
  materialProperties << E, nu, m, n, nKelvin, minTau, timeToDays;

  // instantiate material
  auto material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( LVPLCode,
                                                          materialProperties.data(),
                                                          materialProperties.size(),
                                                          0 ) ) );

  // number of required state vars
  int nStateVars = material->getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  material->assignStateVars( stateVar.data(), nStateVars );

  // declare variables
  Marmot::Vector6d stress;
  Marmot::Vector6d dStrain;
  Marmot::Matrix6d dStressDDStrain;
  double           pNewDT;

  // initialize stress
  stress << 0., 0., 0., 0., 0., 0.;

  // first increment ( load free )
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0;
  time[1] += dT;
  dStrain << 0., 0., 0., 0., 0., 0.;

  // compute material response
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // second increment ( load application )
  dT = 1e-2;
  time[1] += dT;
  double dEps11 = 10e-6;
  dStrain << dEps11, dEps11, 0., 0., 0., 0.;
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  dStrain << 0., 0., 0., 0., 0., 0.;
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // expected stress
  Marmot::Vector6d stressTarget;
  stressTarget << 2.41296e-05, 2.41296e-05, 9.65182e-06, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  Marmot::Testing::throwExceptionOnFailure( Marmot::Testing::checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                                            "Stress computation failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  testLinearViscoelasticPowerLaw();
  return 0;
}
