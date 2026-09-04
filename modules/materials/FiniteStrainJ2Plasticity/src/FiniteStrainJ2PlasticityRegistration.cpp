#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool FiniteStrainJ2PlasticityRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
    FiniteStrainJ2Plasticity >( "FINITESTRAINJ2PLASTICITY" );

  // Register the SUBSTEPPED J2 model
  // This allows you to use "FINITESTRAINJ2PLASTICITY_SUBSTEPPED" in your input file.
  // The properties array must start with [nSubsteps, K, G, fy, ...]
  const static bool FiniteStrainJ2PlasticitySubsteppedRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
    MarmotMaterialFiniteStrainSubstepped< FiniteStrainJ2Plasticity > >( "FINITESTRAINJ2PLASTICITY_SUBSTEPPED" );

} // namespace Marmot::Materials::Registration
