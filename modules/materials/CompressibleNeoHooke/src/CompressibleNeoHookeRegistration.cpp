#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool
      CompressibleNeoHookeRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial< CompressibleNeoHooke >(
        "COMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
