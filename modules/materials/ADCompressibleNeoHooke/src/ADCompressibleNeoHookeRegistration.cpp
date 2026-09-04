#include "Marmot/ADCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool
    ADCompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< ADCompressibleNeoHooke >(
      "ADCOMPRESSIBLENEOHOOKE" );

} // namespace Marmot::Materials::Registration
