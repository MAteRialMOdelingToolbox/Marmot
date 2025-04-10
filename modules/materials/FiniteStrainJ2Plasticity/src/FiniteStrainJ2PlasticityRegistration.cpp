#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

namespace Registration {

constexpr int FiniteStrainJ2PlasticityCode = 11930000 + 16;

using namespace MarmotLibrary;

const static bool FiniteStrainJ2PlasticityRegistered =
    MarmotMaterialFactory::registerMaterial(
        FiniteStrainJ2PlasticityCode, "FINITESTRAINJ2PLASTICITY",
        makeDefaultMarmotMaterialFactoryFunction<
            class FiniteStrainJ2Plasticity>());
} // namespace Registration
} // namespace Marmot::Materials
