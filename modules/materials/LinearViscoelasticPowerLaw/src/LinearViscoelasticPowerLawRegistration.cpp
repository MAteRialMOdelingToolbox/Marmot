#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;


  const static bool LinearViscoelasticPowerLawisRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial<
    LinearViscoelasticPowerLaw >( "LINEARVISCOELASTICPOWERLAW" );

} // namespace
