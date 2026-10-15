#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool ADVonMisesIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< ADVonMises >(
      "ADVONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
