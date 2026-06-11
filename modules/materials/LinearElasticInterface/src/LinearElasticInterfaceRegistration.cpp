#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticInterfaceIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::registerMaterial<
      LinearElasticInterface >( "LINEARELASTICINTERFACE" );

  } // namespace Registration
} // namespace Marmot::Materials
