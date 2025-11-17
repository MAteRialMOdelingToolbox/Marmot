#include "Marmot/Marmot.h"
#include "Marmot/ADVonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool ADVonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial<ADVonMises>(  "ADVONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
