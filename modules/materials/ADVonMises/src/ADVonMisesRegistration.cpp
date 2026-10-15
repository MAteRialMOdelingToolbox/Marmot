#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    ADVonMisesIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< ADVonMises >(
      "ADVONMISES" );

} // namespace
