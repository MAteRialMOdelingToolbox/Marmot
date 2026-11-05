#include "Marmot/LinearViscoElasticWiechert.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool LinearViscoElasticWiechertIsRegistered = Marmot::Registry::MaterialHypoElasticFactory::
    registerMaterial< LinearViscoElasticWiechert >( "LINEARVISCOELASTICWIECHERT" );

} // namespace
