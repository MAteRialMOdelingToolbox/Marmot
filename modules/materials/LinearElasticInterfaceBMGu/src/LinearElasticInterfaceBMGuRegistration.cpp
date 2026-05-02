#include "Marmot/LinearElasticInterfaceBMGu.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial( "LINEARELASTICINTERFACEBMGU",
                        makeDefaultMarmotMaterialHypoElasticInterfaceFactoryFunction<
                          class LinearElasticInterfaceBMGu >() );

  } // namespace Registration
} // namespace Marmot::Materials
