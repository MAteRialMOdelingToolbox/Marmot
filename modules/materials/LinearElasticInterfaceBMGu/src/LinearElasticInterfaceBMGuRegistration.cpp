#include "Marmot/LinearElasticInterfaceBMGu.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticInterfaceBMGuIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial< LinearElasticInterfaceBMGu >( "LINEARELASTICINTERFACEBMGU" );

  } // namespace Registration

} // namespace Marmot::Materials
