#include "Marmot/B4.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

void testB4()
{
  // material properties
  const int b4Code = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( "B4" );
  // material properties
  Eigen::Vector< double, 22 > materialProperties;
  // elastic parameters
  double nu = 0.2;
  double q1 = 20.6;
  // basic creep parameters
  double q2           = 91;
  double q3           = 4.80;
  double q4           = 5.9;
  double n            = 0.1;
  double m            = 0.5;
  int    nKelvinBasic = 12;
  double minTauBasic  = 1e-5;
  // autogenous shrinkage parameters
  double ultimateAutogenousShrinkageStrain = 0.; // -0.0001
  double autogenousShrinkageHalfTime       = 3.;
  double alpha                             = 1.45;
  double rt                                = -4.5;
  // drying shrinkage parameters
  double ultimateDryingShrinkageStrain = -0.0015;
  double dryingShrinkageHalfTime       = 90.;
  double dryingStart                   = 7.;
  double hEnv                          = 1.;
  // drying creep paramters
  double q5            = 400;
  int    nKelvinDrying = 11;
  double minTauDrying  = 2e-4;
  // additional parameters
  double castTime   = -100;
  double timeToDays = 1;
  materialProperties << nu, q1, q2, q3, q4, n, m, nKelvinBasic, minTauBasic, ultimateAutogenousShrinkageStrain,
    autogenousShrinkageHalfTime, alpha, rt, ultimateDryingShrinkageStrain, dryingShrinkageHalfTime, dryingStart, hEnv,
    q5, nKelvinDrying, minTauDrying, castTime, timeToDays;

  // instantiate material
  auto material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( b4Code,
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

  // first increment ( zero strain )
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0;
  time[1] += dT;
  dStrain << 0., 0., 0., 0., 0., 0.;

  // compute material response
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // second increment ( strain application )
  dT = 0.01;
  time[1] += dT;
  dStrain << 10e-6, 10e-6, 0., 0., 0., 0.;
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  dStrain << 0., 0., 0., 0., 0., 0.;
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // expected stress
  Marmot::Vector6d stressTarget;
  stressTarget << 0.377092, 0.377092, 0.150837, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  Marmot::Testing::throwExceptionOnFailure( Marmot::Testing::checkIfEqual< double >( stress, stressTarget, 1e-6 ),
                                            "Stress computation failed for transverse isotropic material in " +
                                              std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  testB4();
  return 0;
}