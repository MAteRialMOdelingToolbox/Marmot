#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool FiniteStrainOrthotropicBiotViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainOrthotropicBiotViscoelasticity >( "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY" );

    const static bool
      FiniteStrainOrthotropicBiotViscoelasticitySubsteppedRegistered = MarmotMaterialFiniteStrainFactory::
        registerMaterial< MarmotMaterialFiniteStrainSubstepped< FiniteStrainOrthotropicBiotViscoelasticity > >(
          "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED" );

  } // namespace Registration
} // namespace Marmot::Materials
