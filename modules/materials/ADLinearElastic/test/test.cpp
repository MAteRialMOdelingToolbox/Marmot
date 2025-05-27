#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot::Testing;

void testADLinearElastic()
{
  // material properties
  const int materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( "ADLINEARELASTIC" );
  Eigen::Vector< double, 2 > materialProperties;
  materialProperties << 210000., 0.3;

  // instantiate material
  auto material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( materialCode,
                                                          materialProperties.data(),
                                                          materialProperties.size(),
                                                          0 ) ) );

  // set increment
  Marmot::Vector6d stress;
  Marmot::Vector6d dStrain;
  Marmot::Matrix6d dStressDDStrain;
  const double     timeOld = 0.0;
  const double     dT      = 1.0;
  double           pNewDT;

  // initialize stress
  stress.setZero();

  // set strain increment
  dStrain << 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825;

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

  // set strain increment
  dStrain << 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825;

  // compute material response
  material->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  std::cout << std::setprecision( 12 ) << stress.transpose() << std::endl;

  // reference solution
  const Marmot::Vector6d stressTarget =
    { 1627.90061538, 416.523692308, -864.896307692, 11.0128846154, 129.673384615, 462.666346154 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-8 ),
                           "comparison with reference solution failed" );
}

int main()
{
  testADLinearElastic();

  return 0;
}
