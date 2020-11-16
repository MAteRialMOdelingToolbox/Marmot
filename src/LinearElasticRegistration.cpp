#include "LinearElastic.h"
#include "MarmotMaterialRegistrationHelper.h"

namespace LinearElasticRegistration {

    using namespace userLibrary;

    const static bool isRegistered = MarmotMaterialFactory::registerMaterial( MaterialCode::LinearElastic,
                                                                              "LINEARELASTIC",
                                                                              makeDefaultMarmotMaterialFactoryFunction<
                                                                                  class LinearElastic >() );

} // namespace LinearElasticRegistration
