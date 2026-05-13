#include "Marmot/LinearElasticInterfaceBMGu.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticInterfaceBMGuIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::
      registerMaterial< LinearElasticInterfaceBMGu >( "LINEARELASTICINTERFACEBMGU" );

  } // namespace Registration

} // namespace Marmot::Materials
