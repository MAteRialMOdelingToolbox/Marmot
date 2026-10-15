#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;


    const static bool LinearElasticIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< LinearElastic >(
      "LINEARELASTIC" );

    // Co-rotational Hughes-Winget wrapper, usable wherever a Marmot::Materials::MaterialFiniteStrain is expected
    // (finite-strain elements, meshfree particles and material points).
    const static bool LinearElasticHughesWingetIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial<
      HughesWingetWrapper< LinearElastic > >( "LINEARELASTIC/HUGHES-WINGET" );

  } // namespace
