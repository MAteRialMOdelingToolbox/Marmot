#include "Marmot/LinearElastic.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial<LinearElastic>( "LINEARELASTIC" );


  } // namespace Registration
} // namespace Marmot::Materials
