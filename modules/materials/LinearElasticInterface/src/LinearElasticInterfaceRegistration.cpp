#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial( "LINEARELASTICINTERFACE",
                        makeDefaultMarmotMaterialHypoElasticInterfaceFactoryFunction<
                          class LinearElasticInterface >() );

  } // namespace Registration
} // namespace Marmot::Materials
