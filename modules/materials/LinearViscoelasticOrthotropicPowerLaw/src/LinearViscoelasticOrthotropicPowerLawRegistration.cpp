#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/Marmot.h"

namespace Marmot::Materials::Registration {

  using namespace MarmotLibrary;

  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = MarmotMaterialFactory::
    registerMaterial<LinearViscoelasticOrthotropicPowerLaw>( "LINEARVISCOELASTICORTHOTROPICPOWERLAW" );

} // namespace Marmot::Materials::Registration
