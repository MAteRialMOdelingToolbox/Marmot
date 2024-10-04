#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int LinearElasticCode = 1;

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( LinearElasticCode,
                        "LINEARELASTIC",
                        makeDefaultMarmotMaterialFactoryFunction< class LinearElastic >() );

  } // namespace Registration
} // namespace Marmot::Materials
