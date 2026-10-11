#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< LinearElastic >(
      "LINEARELASTIC" );

    // Co-rotational Hughes-Winget wrapper, usable wherever a MarmotMaterialFiniteStrain is expected
    // (finite-strain elements, meshfree particles and material points).
    const static bool LinearElasticHughesWingetIsRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      HughesWingetWrapper< LinearElastic > >( "LINEARELASTIC/HUGHES-WINGET" );

  } // namespace Registration
} // namespace Marmot::Materials
