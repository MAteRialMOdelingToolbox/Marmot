#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace {
  using namespace Marmot::Materials;


    const static bool FiniteStrainOrthotropicBiotViscoelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainOrthotropicBiotViscoelasticity >( "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY" );

    const static bool
      FiniteStrainOrthotropicBiotViscoelasticitySubsteppedRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
        registerMaterial< Marmot::Materials::MaterialFiniteStrainSubstepped< FiniteStrainOrthotropicBiotViscoelasticity > >(
          "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED" );

  } // namespace
