#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int vonMisesCode = 2;

    using namespace MarmotLibrary;

    const static bool VonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( vonMisesCode, "VONMISES", makeDefaultMarmotMaterialFactoryFunction< class VonMisesModel >() );

  } // namespace Registration

} // namespace Marmot::Materials
