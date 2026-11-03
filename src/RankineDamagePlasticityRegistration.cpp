#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/RankineDamagePlasticity.h"

namespace {
  using namespace Marmot::Materials;

  const static bool RankineDamagePlasticityIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::
    registerMaterial< RankineDamagePlasticity >( "RANKINEDAMAGEPLASTICITY" );

  // Co-rotational Hughes-Winget wrapper, usable wherever a Marmot::Materials::MaterialFiniteStrain is expected
  const static bool RankineDamagePlasticityHughesWingetIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
    registerMaterial< HughesWingetWrapper< RankineDamagePlasticity > >( "RANKINEDAMAGEPLASTICITY/HUGHES-WINGET" );

} // namespace
