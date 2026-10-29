#include "Marmot/GradientEnhancedFiniteStrainDruckerPrager.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool GradientEnhancedFiniteStrainDruckerPragerRegistered = Marmot::Registry::
    MaterialGradientEnhancedFiniteStrainFactory::registerMaterial< GradientEnhancedFiniteStrainDruckerPrager >(
      "GRADIENTENHANCEDFINITESTRAINDRUCKERPRAGER" );

} // namespace
