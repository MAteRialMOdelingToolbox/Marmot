#include "Marmot/Marmot.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool VonMisesIsRegistered = MarmotMaterialFactory::registerMaterial< VonMisesModel >( "VONMISES" );

  } // namespace Registration

} // namespace Marmot::Materials
