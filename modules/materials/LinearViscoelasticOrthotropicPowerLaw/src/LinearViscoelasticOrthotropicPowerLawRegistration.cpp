#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials::Registration {

  constexpr int LinearViscoelasticOrthotropicPowerLawCode = 1193000 + 21;

  using namespace MarmotLibrary;

  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = MarmotMaterialFactory::
    registerMaterial( LinearViscoelasticOrthotropicPowerLawCode,
                      "LINEARVISCOELASTICORTHOTROPICPOWERLAW",
                      makeDefaultMarmotMaterialFactoryFunction< class LinearViscoelasticOrthotropicPowerLaw >() );

} // namespace Marmot::Materials::Registration
