#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/RankineGradientSofteningPlasticity.h"

namespace {
  using namespace Marmot::Materials;

  const static bool RankineGradientSofteningPlasticityIsRegistered = Marmot::Registry::
    MaterialGradientEnhancedHypoElasticFactory< 1 >::registerMaterial< RankineGradientSofteningPlasticity >(
      "RANKINEGRADIENTSOFTENINGPLASTICITY" );

} // namespace
