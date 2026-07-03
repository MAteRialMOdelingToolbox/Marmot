#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrainFactory.h"

namespace Marmot::Materials {
  const bool isFiniteStrainGradientVonMisesRegistered = MarmotLibrary::
    MarmotMaterialGradientPlasticityFiniteStrainFactory< 1 >::registerMaterial< FiniteStrainGradientVonMises >(
      // Registered in UPPERCASE like GRADIENTVONMISES: the EdelweissMeshfree particle
      // wrapper upper-cases material names before the factory lookup.
      "FINITESTRAINGRADIENTVONMISES" );
}