#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool LinearViscoElasticWiechertIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial<
    LinearViscoElasticWiechert >( "LINEARVISCOELASTICWIECHERT" );

} // namespace Marmot::Materials::Registration
