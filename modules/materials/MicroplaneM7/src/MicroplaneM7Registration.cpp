#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MicroplaneM7.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool MicroplaneM7IsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< MicroplaneM7 >(
      "MICROPLANEM7" );

  } // namespace Registration
} // namespace Marmot::Materials
