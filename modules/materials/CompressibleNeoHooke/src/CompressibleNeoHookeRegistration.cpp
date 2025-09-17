#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int CompressibleNeoHookeCode = 11930000 + 12;

    using namespace MarmotLibrary;

    const static bool CompressibleNeoHookeRegistered = MarmotMaterialFactory::
      registerMaterial( CompressibleNeoHookeCode,
                        "COMPRESSIBLENEOHOOKE",
                        makeDefaultMarmotMaterialFactoryFunction< class CompressibleNeoHooke >() );

  } // namespace Registration
} // namespace Marmot::Materials
