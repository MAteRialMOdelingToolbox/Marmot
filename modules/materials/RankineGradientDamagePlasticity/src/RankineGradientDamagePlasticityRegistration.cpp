#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/RankineGradientDamagePlasticity.h"

namespace {
  using namespace Marmot::Materials;

  const static bool RankineGradientDamagePlasticityIsRegistered = Marmot::Registry::
    MaterialGradientEnhancedHypoElasticFactory< 1 >::registerMaterial< RankineGradientDamagePlasticity >(
      "RANKINEGRADIENTDAMAGEPLASTICITY" );

} // namespace
