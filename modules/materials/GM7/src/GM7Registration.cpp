#include "Marmot/GM7.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotMaterialGradientEnhancedCorotational.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    GM7IsRegistered = Marmot::Registry::MaterialGradientEnhancedHypoElasticFactory< 6 >::registerMaterial< GM7 >(
      "GM7" );

  // GM7 recomputes its stress from its stored strains, so it is made finite-strain in the frame of the polar rotation,
  // not by Hughes-Winget, which hands over a rotated stress that GM7 would overwrite.
  const static bool GM7CorotationalIsRegistered = Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<
    6 >::registerMaterial< GradientEnhancedCorotationalWrapper< GM7 > >( "GM7/COROTATIONAL" );

} // namespace
