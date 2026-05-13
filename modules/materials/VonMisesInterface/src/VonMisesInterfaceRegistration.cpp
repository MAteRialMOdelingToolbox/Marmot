#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/VonMisesInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      VonMisesInterfaceIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::registerMaterial< VonMisesInterface >(
        "VONMISESINTERFACE" );

  } // namespace Registration

} // namespace Marmot::Materials
