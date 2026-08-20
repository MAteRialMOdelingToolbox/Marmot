#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/Ogden.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool OgdenRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< Ogden >( "OGDEN" );

  } // namespace Registration
} // namespace Marmot::Materials
