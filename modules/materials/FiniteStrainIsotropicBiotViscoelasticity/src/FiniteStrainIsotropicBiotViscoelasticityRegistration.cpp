#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool FiniteStrainIsotropicBiotViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainIsotropicBiotViscoelasticity >( "FINITESTRAINISOTROPICBIOTVISCOELASTICITY" );

  } // namespace Registration
} // namespace Marmot::Materials
