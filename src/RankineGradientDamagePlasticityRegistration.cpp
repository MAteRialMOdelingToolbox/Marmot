#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/RankineGradientDamagePlasticity.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      RankineGradientDamagePlasticityIsRegistered = MarmotMaterialGeneralGradientEnhancedHypoElasticFactory<
        1 >::registerMaterial< RankineGradientDamagePlasticity >( "RANKINEGRADIENTDAMAGEPLASTICITY" );

  } // namespace Registration

} // namespace Marmot::Materials
