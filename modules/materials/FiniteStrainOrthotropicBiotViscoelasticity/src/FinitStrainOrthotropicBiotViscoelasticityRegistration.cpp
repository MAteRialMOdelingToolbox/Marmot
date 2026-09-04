#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool FiniteStrainOrthotropicBiotViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< FiniteStrainOrthotropicBiotViscoelasticity >( "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY" );

  const static bool FiniteStrainOrthotropicBiotViscoelasticitySubsteppedRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< MarmotMaterialFiniteStrainSubstepped< FiniteStrainOrthotropicBiotViscoelasticity > >(
      "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED" );

} // namespace Marmot::Materials::Registration
