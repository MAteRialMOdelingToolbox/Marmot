#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {


  const static bool LinearViscoElasticWiechertIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial<
    LinearViscoElasticWiechert >( "LINEARVISCOELASTICWIECHERT" );

} // namespace Marmot::Materials::Registration
