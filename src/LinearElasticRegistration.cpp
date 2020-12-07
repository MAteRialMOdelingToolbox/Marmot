#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace LinearElasticRegistration {

    using namespace MarmotLibrary;

    const static bool isRegistered = MarmotMaterialFactory::registerMaterial( MaterialCode::LinearElastic,
                                                                              "LINEARELASTIC",
                                                                              makeDefaultMarmotMaterialFactoryFunction<
                                                                                  class LinearElastic >() );

} // namespace LinearElasticRegistration
