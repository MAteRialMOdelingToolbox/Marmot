#include "Marmot/B4.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool
    B4isRegistered = MarmotMaterialFactory::registerMaterial( MaterialCode::B4,
                                                              "B4",
                                                              makeDefaultMarmotMaterialFactoryFunction< class B4 >() );

} // namespace Marmot::Materials::Registration
