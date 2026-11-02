#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = Marmot::Registry::MaterialHypoElasticFactory::
    registerMaterial< LinearViscoelasticOrthotropicPowerLaw >( "LINEARVISCOELASTICORTHOTROPICPOWERLAW" );

} // namespace
