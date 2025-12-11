#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      ADLinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< ADLinearElastic >(
        "ADLINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
