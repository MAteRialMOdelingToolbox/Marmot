#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool
      GradientEnhancedFiniteStrainDruckerPragerRegistered = Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory::
        registerMaterial< GradientEnhancedFiniteStrainDruckerPrager >( "GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER" );

  } // namespace Registration
} // namespace Marmot::Materials
