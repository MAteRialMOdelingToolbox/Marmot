#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool LinearViscoelasticPowerLawisRegistered = MarmotMaterialHypoElasticFactory::registerMaterial<
    LinearViscoelasticPowerLaw >( "LINEARVISCOELASTICPOWERLAW" );

} // namespace Marmot::Materials::Registration
