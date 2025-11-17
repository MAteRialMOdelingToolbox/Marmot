#include "Marmot/Marmot.h"
#include "Marmot/B4.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool
    B4isRegistered = MarmotMaterialFactory::registerMaterial<B4>( "B4" );

} // namespace Marmot::Materials::Registration
