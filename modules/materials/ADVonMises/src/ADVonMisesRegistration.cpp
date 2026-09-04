#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool ADVonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< ADVonMises >(
    "ADVONMISES" );

} // namespace Marmot::Materials::Registration
