#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {


  const static bool LinearViscoelasticPowerLawisRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial<
    LinearViscoelasticPowerLaw >( "LINEARVISCOELASTICPOWERLAW" );

} // namespace Marmot::Materials::Registration
