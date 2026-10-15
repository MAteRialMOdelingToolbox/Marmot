#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool CompressibleFiniteStrainLinearViscoelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< CompressibleFiniteStrainLinearViscoelasticity >(
        "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY" );

    const static bool
      CompressibleFiniteStrainLinearViscoelasticitySubsteppedRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
        registerMaterial< Marmot::Materials::MaterialFiniteStrainSubstepped< CompressibleFiniteStrainLinearViscoelasticity > >(
          "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED" );
  } // namespace Registration
} // namespace Marmot::Materials
