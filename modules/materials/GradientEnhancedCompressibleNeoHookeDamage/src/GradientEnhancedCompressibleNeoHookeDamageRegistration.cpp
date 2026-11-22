#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool GradientEnhancedCompressibleNeoHookeDamageRegistered = Marmot::Registry::
    MaterialGradientEnhancedFiniteStrainFactory< 1 >::registerMaterial< GradientEnhancedCompressibleNeoHookeDamage >(
      "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE" );

} // namespace
