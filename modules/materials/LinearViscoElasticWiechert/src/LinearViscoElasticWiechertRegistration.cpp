#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticWiechertIsRegistered =
      MarmotMaterialHypoElasticFactory::registerMaterial< LinearViscoElasticWiechert >(
        "LINEARVISCOELASTICWIECHERT" );

  } // namespace Registration

} // namespace Marmot::Materials
