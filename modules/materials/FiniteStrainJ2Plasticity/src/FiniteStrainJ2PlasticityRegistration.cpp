#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool FiniteStrainJ2PlasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial<
      FiniteStrainJ2Plasticity >( "FINITESTRAINJ2PLASTICITY" );

    // Register the SUBSTEPPED J2 model
    // This allows you to use "FINITESTRAINJ2PLASTICITY_SUBSTEPPED" in your input file.
    // The properties array must start with [nSubsteps, K, G, fy, ...]
    const static bool FiniteStrainJ2PlasticitySubsteppedRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< Marmot::Materials::MaterialFiniteStrainSubstepped< FiniteStrainJ2Plasticity > >(
        "FINITESTRAINJ2PLASTICITY_SUBSTEPPED" );

  } // namespace Registration
} // namespace Marmot::Materials
