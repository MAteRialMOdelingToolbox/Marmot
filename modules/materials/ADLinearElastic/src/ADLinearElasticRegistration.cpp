#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool ADLinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< ADLinearElastic >(
    "ADLINEARELASTIC" );

} // namespace Marmot::Materials::Registration
