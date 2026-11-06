#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/RankineSofteningPlasticity.h"

namespace {
  using namespace Marmot::Materials;

  const static bool RankineSofteningPlasticityIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::
    registerMaterial< RankineSofteningPlasticity >( "RANKINESOFTENINGPLASTICITY" );

  // Co-rotational Hughes-Winget wrapper, usable wherever a Marmot::Materials::MaterialFiniteStrain is expected
  // (finite-strain elements, meshfree particles and material points).
  const static bool RankineSofteningPlasticityHughesWingetIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
    registerMaterial< HughesWingetWrapper< RankineSofteningPlasticity > >( "RANKINESOFTENINGPLASTICITY/HUGHES-WINGET" );

} // namespace
