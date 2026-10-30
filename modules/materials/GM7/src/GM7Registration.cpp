#include "Marmot/GM7.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      GM7IsRegistered = MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::registerMaterial< GM7 >( "GM7" );

  } // namespace Registration
} // namespace Marmot::Materials
