#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/VonMisesInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      VonMisesInterfaceIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::registerMaterial< VonMisesInterface >(
        "VONMISESINTERFACE" );

  } // namespace Registration

} // namespace Marmot::Materials
