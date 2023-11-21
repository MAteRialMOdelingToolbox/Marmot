#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials::Registration {
  
  constexpr int LinearViscoelasticPowerLawCode = 1193000 + 20;

  using namespace MarmotLibrary;

  const static bool
    LinearViscoelasticPowerLawisRegistered = MarmotMaterialFactory::registerMaterial( LinearViscoelasticPowerLawCode,
                                                              "LINEARVISCOELASTICPOWERLAW",
                                                              makeDefaultMarmotMaterialFactoryFunction< class LinearViscoelasticPowerLaw >() );

} // namespace Marmot::Materials::Registration
