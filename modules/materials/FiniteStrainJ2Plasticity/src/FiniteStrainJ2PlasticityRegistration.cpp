#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool FiniteStrainJ2PlasticityRegistered = MarmotMaterialFactory::
      registerMaterial<FiniteStrainJ2Plasticity >( "FINITESTRAINJ2PLASTICITY");
  } // namespace Registration
} // namespace Marmot::Materials
