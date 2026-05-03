#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticInterfaceIsRegistered = MarmotMaterialHypoElasticInterfaceFactory::
      registerMaterial< LinearViscoElasticInterface >( "LINEARVISCOELASTICINTERFACE" );

  } // namespace Registration

} // namespace Marmot::Materials
