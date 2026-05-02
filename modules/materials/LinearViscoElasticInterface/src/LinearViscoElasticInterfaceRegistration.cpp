#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial( "LINEARVISCOELASTICINTERFACE",
                        makeDefaultMarmotMaterialHypoElasticInterfaceFactoryFunction<
                          class LinearViscoElasticInterface >() );

  } // namespace Registration
} // namespace Marmot::Materials
