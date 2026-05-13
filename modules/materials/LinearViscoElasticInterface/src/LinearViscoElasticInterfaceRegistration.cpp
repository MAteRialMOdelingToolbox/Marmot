#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticInterfaceIsRegistered = MarmotInterfaceMaterialHypoElasticFactory::
      registerMaterial< LinearViscoElasticInterface >( "LINEARVISCOELASTICINTERFACE" );

  } // namespace Registration

} // namespace Marmot::Materials
