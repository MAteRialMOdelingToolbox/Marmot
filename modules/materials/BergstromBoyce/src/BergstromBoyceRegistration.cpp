#include "Marmot/BergstromBoyce.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    BergstromBoyceRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial< BergstromBoyce >(
      "BERGSTROMBOYCE" );
} // namespace
