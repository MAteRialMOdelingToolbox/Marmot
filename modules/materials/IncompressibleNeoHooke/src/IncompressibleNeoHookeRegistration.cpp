#include "Marmot/IncompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      IncompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< IncompressibleNeoHooke >(
        "INCOMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
