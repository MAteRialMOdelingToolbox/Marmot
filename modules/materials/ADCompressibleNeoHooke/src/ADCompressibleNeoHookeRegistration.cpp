#include "Marmot/ADCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;


    const static bool
      ADCompressibleNeoHookeRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial< ADCompressibleNeoHooke >(
        "ADCOMPRESSIBLENEOHOOKE" );

  } // namespace
