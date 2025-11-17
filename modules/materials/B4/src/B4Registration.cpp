#include "Marmot/B4.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool B4isRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< B4 >( "B4" );

} // namespace Marmot::Materials::Registration
