#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool FiniteStrainJ2PlasticityRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      FiniteStrainJ2Plasticity >( "FINITESTRAINJ2PLASTICITY" );
  } // namespace Registration
} // namespace Marmot::Materials
