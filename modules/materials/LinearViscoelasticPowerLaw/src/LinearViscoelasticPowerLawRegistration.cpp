#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool LinearViscoelasticPowerLawisRegistered = MarmotMaterialFactory::
    registerMaterial<LinearViscoelasticPowerLaw>( 
                      "LINEARVISCOELASTICPOWERLAW"
                       );

} // namespace Marmot::Materials::Registration
