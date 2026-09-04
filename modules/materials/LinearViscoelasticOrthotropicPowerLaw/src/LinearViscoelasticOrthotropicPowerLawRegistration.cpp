#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = MarmotMaterialHypoElasticFactory::
    registerMaterial< LinearViscoelasticOrthotropicPowerLaw >( "LINEARVISCOELASTICORTHOTROPICPOWERLAW" );

} // namespace Marmot::Materials::Registration
