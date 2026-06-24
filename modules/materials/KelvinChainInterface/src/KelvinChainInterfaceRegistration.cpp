#include "Marmot/KelvinChainInterface.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool KelvinChainInterfaceIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::registerMaterial<
    KelvinChainInterface >( "KELVINCHAININTERFACE" );

} // namespace Marmot::Materials::Registration
