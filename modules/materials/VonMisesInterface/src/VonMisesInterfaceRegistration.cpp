#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/VonMisesInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool VonMisesIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial( "VONMISESINTERFACE",
                        makeDefaultMarmotMaterialHypoElasticInterfaceFactoryFunction< class VonMisesInterface >() );

  } // namespace Registration
} // namespace Marmot::Materials
