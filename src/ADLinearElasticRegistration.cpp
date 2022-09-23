#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADLinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( MaterialCode::ADLinearElastic,
                        "ADLINEARELASTIC",
                        makeDefaultMarmotMaterialFactoryFunction< class ADLinearElastic >() );

  } // namespace Registration
} // namespace Marmot::Materials
