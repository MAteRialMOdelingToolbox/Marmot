#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialHughesWinget.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/VonMises.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    VonMisesIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< VonMisesModel >(
      "VONMISES" );

  // Co-rotational Hughes-Winget wrapper, usable wherever a Marmot::Materials::MaterialFiniteStrain is expected
  // (finite-strain elements, meshfree particles and material points).
  const static bool VonMisesHughesWingetIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::registerMaterial<
    HughesWingetWrapper< VonMisesModel > >( "VONMISES/HUGHES-WINGET" );

  // Variant recovering the exact d(sigma^(n+1))/d(sigma_rot) operator; worth the seven material
  // evaluations whenever large incremental rotations coincide with plastic flow.
  const static bool VonMisesHughesWingetExactIsRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
    registerMaterial< HughesWingetWrapper< VonMisesModel, HughesWingetTangent::Exact > >(
      "VONMISES/HUGHES-WINGET/EXACT-TANGENT" );

} // namespace
