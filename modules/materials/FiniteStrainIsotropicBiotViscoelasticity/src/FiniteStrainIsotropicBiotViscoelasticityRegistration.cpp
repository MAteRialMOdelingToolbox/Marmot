#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;


    const static bool FiniteStrainIsotropicBiotViscoelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
      registerMaterial< FiniteStrainIsotropicBiotViscoelasticity >( "FINITESTRAINISOTROPICBIOTVISCOELASTICITY" );

  } // namespace
