#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/WiechertInterfaceMaterial.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool WiechertInterfaceMaterialIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::registerMaterial<
    WiechertInterfaceMaterial >( "WIECHERTINTERFACE" );

} // namespace Marmot::Materials::Registration
