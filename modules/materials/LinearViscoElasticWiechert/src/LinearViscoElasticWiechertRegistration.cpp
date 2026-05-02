#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {
    constexpr int LinearViscoelasticWiechertCode = 1193000 + 22222;

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( LinearViscoelasticWiechertCode,
                        "LINEARVISCOELASTICWIECHERT",
                        makeDefaultMarmotMaterialFactoryFunction< class LinearViscoElasticWiechert >() );

  } // namespace Registration
} // namespace Marmot::Materials
