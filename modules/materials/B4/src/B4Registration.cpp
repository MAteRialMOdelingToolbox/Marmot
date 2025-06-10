#include "Marmot/B4.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials::Registration {

  constexpr int B4Code = 1193001;

  using namespace MarmotLibrary;

  const static bool
    B4isRegistered = MarmotMaterialFactory::registerMaterial( B4Code,
                                                              "B4",
                                                              makeDefaultMarmotMaterialFactoryFunction< class B4 >() );

} // namespace Marmot::Materials::Registration
