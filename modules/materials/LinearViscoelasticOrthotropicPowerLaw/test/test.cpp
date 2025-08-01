#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"

#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

void testLinearViscoelasticOrthotropicPowerLaw()
{

  double E  = 2e5;
  double nu = 0.2;
  double G  = E / ( 2 * ( 1 + nu ) );
  // material properties
  double materialProperties[22];
  // Young's moduli
  materialProperties[0] = E;
  materialProperties[1] = E;
  materialProperties[2] = E;
  // Poisson's ratios
  materialProperties[3] = nu;
  materialProperties[4] = nu;
  materialProperties[5] = nu;
  // shear moduli
  materialProperties[6] = G;
  materialProperties[7] = G;
  materialProperties[8] = G;

  // viscoelastic parameters
  materialProperties[9]  = 0.5;
  materialProperties[10] = 0.1;
  materialProperties[11] = 2;
  materialProperties[12] = 10;
  materialProperties[13] = 0.0001;
  materialProperties[14] = 3.1622776601683795;
  materialProperties[15] = 1.;
  // coordinate system
  materialProperties[16] = 1.;
  materialProperties[17] = 0.5;
  materialProperties[18] = 0.;
  materialProperties[19] = -0.5;
  materialProperties[20] = 1;
  materialProperties[21] = 1;

  // instantiate material
  auto material = Marmot::Materials::LinearViscoelasticOrthotropicPowerLaw( &materialProperties[0], 20, 1 );

  // number of required state vars
  int nStateVars = material.getNumberOfRequiredStateVars();

  // initialize state vars
  Eigen::VectorXd stateVar( nStateVars );
  stateVar.setZero();
  material.assignStateVars( stateVar.data(), nStateVars );

  // declare variables
  Marmot::Vector6d stress;
  Marmot::Vector6d dStrain;
  Marmot::Matrix6d dStressDDStrain;
  double           pNewDT;

  // initialize stress
  stress.setZero();

  // first increment ( load free )
  Eigen::VectorXd time( 2 );
  time.setZero();
  double dT = 28.0;
  time[1] += dT;
  dStrain.setZero();
  // compute material response
  material.computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // second increment ( load application )
  dT = 1e-2;
  time[1] += dT;
  double dEps11 = 10e-6;
  // dStrain.setConstant( dEps11 );
  dStrain << dEps11, dEps11, 0., 0., 0., 0.;
  material.computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // third increment ( constant strain, relaxation )
  dT = 100.;
  time[1] += dT;
  dStrain << 0., 0., 0., 0., 0., 0.;
  material.computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), time.data(), dT, pNewDT );

  // expected stress
  Marmot::Vector6d stressTarget;
  stressTarget << 2.45582e-05, 2.45582e-05, 9.82329e-06, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  Marmot::Testing::throwExceptionOnFailure( Marmot::Testing::checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                                            "Stress computation failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  testLinearViscoelasticOrthotropicPowerLaw();
  return 0;
}
