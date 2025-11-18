#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADVonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< ADVonMises >(
      "ADVONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
