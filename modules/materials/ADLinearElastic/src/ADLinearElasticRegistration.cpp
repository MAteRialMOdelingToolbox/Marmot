#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include <cstddef>
#include <string>

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    constexpr int base                = 11930000;
    constexpr int code                = 1;
    constexpr int ADLinearElasticCode = base + code;

    const static bool ADLinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( ADLinearElasticCode,
                        "ADLINEARELASTIC",
                        makeDefaultMarmotMaterialFactoryFunction< class ADLinearElastic >() );

  } // namespace Registration
} // namespace Marmot::Materials
