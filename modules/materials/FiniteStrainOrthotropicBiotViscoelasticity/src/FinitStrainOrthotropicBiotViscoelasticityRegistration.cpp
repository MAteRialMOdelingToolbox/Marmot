#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool FiniteStrainOrthotropicBiotViscoelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainOrthotropicBiotViscoelasticity >( "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY" );

    const static bool
      FiniteStrainOrthotropicBiotViscoelasticitySubsteppedRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
        registerMaterial< Marmot::Materials::MaterialFiniteStrainSubstepped< FiniteStrainOrthotropicBiotViscoelasticity > >(
          "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED" );

  } // namespace Registration
} // namespace Marmot::Materials
