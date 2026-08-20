#include "Marmot/IncompressibleMooneyRivlin.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool IncompressibleMooneyRivlinRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial<
      IncompressibleMooneyRivlin >( "INCOMPRESSIBLEMOONEYRIVLIN" );

  } // namespace Registration
} // namespace Marmot::Materials
