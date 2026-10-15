#include "Marmot/ADCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool
      ADCompressibleNeoHookeRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial< ADCompressibleNeoHooke >(
        "ADCOMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
