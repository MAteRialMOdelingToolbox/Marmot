#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool FiniteStrainIsotropicBiotViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< FiniteStrainIsotropicBiotViscoelasticity >( "FINITESTRAINISOTROPICBIOTVISCOELASTICITY" );

} // namespace Marmot::Materials::Registration
