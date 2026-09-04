#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool CompressibleFiniteStrainLinearViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< CompressibleFiniteStrainLinearViscoelasticity >(
      "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY" );

  const static bool
    CompressibleFiniteStrainLinearViscoelasticitySubsteppedRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< MarmotMaterialFiniteStrainSubstepped< CompressibleFiniteStrainLinearViscoelasticity > >(
        "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED" );
} // namespace Marmot::Materials::Registration
