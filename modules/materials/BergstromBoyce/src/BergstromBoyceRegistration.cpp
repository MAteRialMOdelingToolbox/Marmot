#include "Marmot/BergstromBoyce.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool BergstromBoyceRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      BergstromBoyce >( "BERGSTROMBOYCE" );

  } // namespace Registration
} // namespace Marmot::Materials
