#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

    namespace Registration {

        using namespace MarmotLibrary;

        const static bool LinearElasticIsRegistered = MarmotMaterialFactory::
            registerMaterial( MaterialCode::LinearElastic,
                              "LINEARELASTIC",
                              makeDefaultMarmotMaterialFactoryFunction< class LinearElastic >() );

    } // namespace Registration
} // namespace Marmot::Materials
