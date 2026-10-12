#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool VonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< VonMisesModel >(
      "VONMISES" );

    // Co-rotational Hughes-Winget wrapper, usable wherever a MarmotMaterialFiniteStrain is expected
    // (finite-strain elements, meshfree particles and material points).
    const static bool VonMisesHughesWingetIsRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      HughesWingetWrapper< VonMisesModel > >( "VONMISES/HUGHES-WINGET" );

    // Variant recovering the exact d(sigma^(n+1))/d(sigma_rot) operator; worth the seven material
    // evaluations whenever large incremental rotations coincide with plastic flow.
    const static bool VonMisesHughesWingetExactIsRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      HughesWingetWrapper< VonMisesModel, HughesWingetTangent::Exact > >( "VONMISES/HUGHES-WINGET/EXACT-TANGENT" );

  } // namespace Registration

} // namespace Marmot::Materials
