#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool CompressibleNeoHookeRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial<
    CompressibleNeoHooke >( "COMPRESSIBLENEOHOOKE" );

} // namespace
