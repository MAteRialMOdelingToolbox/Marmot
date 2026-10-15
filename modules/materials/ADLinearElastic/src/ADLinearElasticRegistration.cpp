#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool
      ADLinearElasticIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< ADLinearElastic >(
        "ADLINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
