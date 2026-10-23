#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/RankineDamagePlasticity.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool RankineDamagePlasticityIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial<
      RankineDamagePlasticity >( "RANKINEDAMAGEPLASTICITY" );

    // Co-rotational Hughes-Winget wrapper, usable wherever a MarmotMaterialFiniteStrain is expected
    const static bool RankineDamagePlasticityHughesWingetIsRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< HughesWingetWrapper< RankineDamagePlasticity > >( "RANKINEDAMAGEPLASTICITY/HUGHES-WINGET" );

  } // namespace Registration

} // namespace Marmot::Materials
