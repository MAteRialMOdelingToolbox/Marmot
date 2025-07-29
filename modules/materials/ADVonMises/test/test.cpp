#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot::Testing;

void testADVonMises()
{
  // material properties
  const int vonMisesCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( "ADVONMISES" );
  Eigen::Vector< double, 6 > materialProperties;
  materialProperties << 210000., 0.3, 200., 2100., 20., 20;

  // instantiate material
  auto material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( vonMisesCode,
                                                          materialProperties.data(),
                                                          materialProperties.size(),
                                                          0 ) ) );
  // set state vars
  if ( material->getNumberOfRequiredStateVars() > 1 ) {
    throw std::runtime_error( "Number of required state vars changed!" );
  }

  // set increment
  Marmot::Vector6d stress;
  Marmot::Vector6d dStrain;
  Marmot::Matrix6d dStressDDStrain;
  const double     timeOld = 0.0;
  const double     dT      = 1.0;
  double           pNewDT;

  // initialize stress
  stress.setZero();

  // initialize material state
  double kappa = 0.;
  material->assignStateVars( &kappa, 1 );

  // set strain increment
  dStrain << 0.01, 0, 0., 0., 0.03, 0.;

  // check coordinate invariance
  throwExceptionOnFailure( spinTurbokreisel( material,
                                             stress.data(),
                                             dStressDDStrain.data(),
                                             dStrain.data(),
                                             &timeOld,
                                             dT,
                                             pNewDT,
                                             1e-10,
                                             1e-8 ),
                           "Turbokreisel failed!" );

  // initialize stress
  stress.setZero();

  // initialize material state
  kappa = 0.;
  material->assignStateVars( &kappa, 1 );

  // set strain increment
  dStrain << 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825;

  // compute material response
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // reference solution
  const Marmot::Vector6d stressTarget =
    { 511.262747695, 395.408929527, 272.856322779, 1.0532516407, 12.4017194288, 44.2485420671 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-9 ),
                           "comparison with reference solution failed" );
}

int main()
{
  testADVonMises();

  return 0;
}
