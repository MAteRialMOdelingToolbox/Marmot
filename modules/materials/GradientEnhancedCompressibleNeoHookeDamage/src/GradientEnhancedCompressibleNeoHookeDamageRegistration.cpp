#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool
      GradientEnhancedCompressibleNeoHookeDamageRegistered = Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory::
        registerMaterial< GradientEnhancedCompressibleNeoHookeDamage >( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE" );

  } // namespace Registration
} // namespace Marmot::Materials
