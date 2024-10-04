#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "memory"
#include <array>
#include <iostream>

int main()
{

  const std::array< double, 2 > materialProperties  = { 20000, 0.25 };
  const double                  nMaterialProperties = 2;
  const int                     elLabel             = 1;

  std::unique_ptr< MarmotMaterialHypoElastic > mat( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName(
                                                            "LINEARELASTIC" ),
                                                          materialProperties.data(),
                                                          nMaterialProperties,
                                                          elLabel ) ) );

  Marmot::Vector6d stress  = { 0., 0., 0., 0., 0., 0. };
  Marmot::Vector6d dStrain = { 1e-3, 0., 0., 0., 0., 0. };
  Marmot::Matrix6d dStressDDStrain;
  const double     timeOld = 0.0;
  const double     dT      = 1.0;
  double           pNewDT;

  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  const Marmot::Vector6d stressTarget = { 24., 8., 8., 0., 0., 0. };

  if ( ( stress - stressTarget ).norm() > 1e-10 )
    return 1;

  return 0;
}
