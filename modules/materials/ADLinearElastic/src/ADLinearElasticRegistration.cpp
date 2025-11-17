#include "Marmot/ADLinearElastic.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADLinearElasticIsRegistered = MarmotMaterialFactory::registerMaterial< ADLinearElastic >(
      "ADLINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
