#include "Marmot/GM7.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    GM7IsRegistered = Marmot::Registry::MaterialGradientEnhancedHypoElasticFactory< 6 >::registerMaterial< GM7 >(
      "GM7" );

} // namespace
