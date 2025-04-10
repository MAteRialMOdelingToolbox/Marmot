#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    constexpr int base           = 11930000;
    constexpr int code           = 3;
    constexpr int ADVonMisesCode = base + code;

    const static bool ADVonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( ADVonMisesCode, "ADVONMISES", makeDefaultMarmotMaterialFactoryFunction< class ADVonMises >() );

  } // namespace Registration
} // namespace Marmot::Materials
