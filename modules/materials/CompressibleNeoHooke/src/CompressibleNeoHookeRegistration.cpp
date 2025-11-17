#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool CompressibleNeoHookeRegistered = MarmotMaterialFactory::registerMaterial< CompressibleNeoHooke >(
      "COMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
