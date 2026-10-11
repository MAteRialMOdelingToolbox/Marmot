#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      GradientEnhancedCompressibleNeoHookeDamageRegistered = MarmotMaterialGradientEnhancedFiniteStrainFactory::
        registerMaterial< GradientEnhancedCompressibleNeoHookeDamage >( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE" );

  } // namespace Registration
} // namespace Marmot::Materials
