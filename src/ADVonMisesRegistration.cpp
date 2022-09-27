#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADVonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( MaterialCode::ADVonMises,
                        "ADVONMISES",
                        makeDefaultMarmotMaterialFactoryFunction< class ADVonMises >() );

  } // namespace Registration
} // namespace Marmot::Materials
