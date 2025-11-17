#include "Marmot/ADVonMises.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADVonMisesIsRegistered = MarmotMaterialFactory::registerMaterial< ADVonMises >( "ADVONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
