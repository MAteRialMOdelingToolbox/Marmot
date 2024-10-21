#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "memory"
#include <array>

std::unique_ptr< MarmotMaterialHypoElastic > createMarmotMaterialHypoElastic( std::string   materialName,
                                                                              const double* materialProperties,
                                                                              const int     nMaterialProperties )
{

  const int                                    elLabel = 1;
  std::unique_ptr< MarmotMaterialHypoElastic > mat( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName(
                                                            materialName ),
                                                          materialProperties,
                                                          nMaterialProperties,
                                                          elLabel ) ) );
  return mat;
}

int main()
{
  // initialize material parameters
  const std::array< double, 2 > materialProperties  = { 20000, 0.25 };
  const double                  nMaterialProperties = 2;

  // create the material
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties.data(), nMaterialProperties );

  // set increment
  Marmot::Vector6d stress  = { 0., 0., 0., 0., 0., 0. };
  Marmot::Vector6d dStrain = { 1e-3, 0., 0., 0., 0., 0. };
  Marmot::Matrix6d dStressDDStrain;
  const double     timeOld = 0.0;
  const double     dT      = 1.0;
  double           pNewDT;

  // compute material response
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // set the correct result
  const Marmot::Vector6d stressTarget = { 24., 8., 8., 0., 0., 0. };

  // check the result
  Marmot::Testing::checkIfEqual< double >( stress, stressTarget );

  // TODO: add tests also for shear increments and anisotropic formulations

  return 0;
}
