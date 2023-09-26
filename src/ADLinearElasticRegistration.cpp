#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    constexpr int ADLinearElasticCode = 1193002;

    const static bool ADLinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( ADLinearElasticCode,
                        "ADLINEARELASTIC",
                        makeDefaultMarmotMaterialFactoryFunction< class ADLinearElastic >() );

  } // namespace Registration
} // namespace Marmot::Materials
