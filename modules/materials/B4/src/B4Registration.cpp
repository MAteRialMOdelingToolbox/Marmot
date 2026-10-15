#include "Marmot/B4.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace {
  using namespace Marmot::Materials;


  const static bool B4isRegistered = Marmot::Registry::MaterialHypoElasticFactory::registerMaterial< B4 >( "B4" );

} // namespace
