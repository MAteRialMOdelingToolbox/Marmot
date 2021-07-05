#include "Marmot/B3.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool B3isRegistered = MarmotMaterialFactory::
    registerMaterial( MaterialCode::B3,
                      "B3",
                      makeDefaultMarmotMaterialFactoryFunction< class B3 >() );

} // namespace Marmot::Materials::Registration
