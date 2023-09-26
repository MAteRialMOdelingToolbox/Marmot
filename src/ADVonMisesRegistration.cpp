#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    constexpr int ADVonMisesCode = 1193003;

    const static bool ADVonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( ADVonMisesCode, "ADVONMISES", makeDefaultMarmotMaterialFactoryFunction< class ADVonMises >() );

  } // namespace Registration
} // namespace Marmot::Materials
