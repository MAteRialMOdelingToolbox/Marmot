#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< LinearElastic >(
      "LINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
