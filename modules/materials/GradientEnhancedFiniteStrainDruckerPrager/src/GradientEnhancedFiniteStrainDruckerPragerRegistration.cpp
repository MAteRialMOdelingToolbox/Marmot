#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      GradientEnhancedFiniteStrainDruckerPragerRegistered = MarmotMaterialGradientEnhancedFiniteStrainFactory::
        registerMaterial< GradientEnhancedFiniteStrainDruckerPrager >( "GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER" );

  } // namespace Registration
} // namespace Marmot::Materials
