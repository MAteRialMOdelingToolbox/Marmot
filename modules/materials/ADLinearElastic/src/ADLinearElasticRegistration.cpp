#include "Marmot/Marmot.h"
#include "Marmot/ADLinearElastic.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADLinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial<ADLinearElastic>( "ADLINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
