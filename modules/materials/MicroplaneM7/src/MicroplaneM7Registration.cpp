#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MicroplaneM7.h"

namespace {
  using namespace Marmot::Materials;

  const static bool
    MicroplaneM7IsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< MicroplaneM7 >(
      "MICROPLANEM7" );

} // namespace
