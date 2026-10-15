#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {


  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = Marmot::Registry::MaterialHypoElasticFactory::
    registerMaterial< LinearViscoelasticOrthotropicPowerLaw >( "LINEARVISCOELASTICORTHOTROPICPOWERLAW" );

} // namespace Marmot::Materials::Registration
