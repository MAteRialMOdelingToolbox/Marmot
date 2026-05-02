#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticInterfaceIsRegistered =
      MarmotMaterialHypoElasticInterfaceFactory::registerMaterial< LinearElasticInterface >( "LINEARELASTICINTERFACE" );

  } // namespace Registration

} // namespace Marmot::Materials
