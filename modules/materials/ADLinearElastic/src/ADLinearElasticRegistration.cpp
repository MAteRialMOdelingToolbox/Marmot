#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;


    const static bool
      ADLinearElasticIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< ADLinearElastic >(
        "ADLINEARELASTIC" );

  } // namespace
