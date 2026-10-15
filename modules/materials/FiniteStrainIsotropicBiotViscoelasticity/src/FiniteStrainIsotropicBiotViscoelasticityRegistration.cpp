#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {


    const static bool FiniteStrainIsotropicBiotViscoelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainIsotropicBiotViscoelasticity >( "FINITESTRAINISOTROPICBIOTVISCOELASTICITY" );

  } // namespace Registration
} // namespace Marmot::Materials
