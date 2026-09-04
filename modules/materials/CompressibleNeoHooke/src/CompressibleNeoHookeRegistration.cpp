#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool
    CompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< CompressibleNeoHooke >(
      "COMPRESSIBLENEOHOOKE" );

} // namespace Marmot::Materials::Registration
