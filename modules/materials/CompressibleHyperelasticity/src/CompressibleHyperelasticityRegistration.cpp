#include "Marmot/CompressibleHyperelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool CompressibleHyperelasticityRegistered = Marmot::Registry::MaterialFiniteStrainFactory::
    registerMaterial< CompressibleHyperelasticity >( "COMPRESSIBLEHYPERELASTICITY" );

} // namespace
