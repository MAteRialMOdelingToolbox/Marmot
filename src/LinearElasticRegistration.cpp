#include "LinearElastic.h"
#include "bftMaterialRegistrationHelper.h"

namespace LinearElasticRegistration {

    using namespace userLibrary;

    const static bool
        isRegistered = BftMaterialFactory::registerMaterial( MaterialCode::LinearElastic,
                                                             "LINEARELASTIC",
                                                             makeDefaultBftMaterialFactoryFunction<class LinearElastic>() );

} // namespace LinearElasticRegistration
