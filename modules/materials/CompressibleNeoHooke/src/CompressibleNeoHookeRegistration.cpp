#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      CompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< CompressibleNeoHooke >(
        "COMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
