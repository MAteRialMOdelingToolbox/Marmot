#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

extern "C" __attribute__( ( constructor, used, visibility( "default" ) ) ) void MarmotRegisterLinearElasticInterface()
{
  static bool registered = MarmotLibrary::MarmotInterfaceMaterialHypoElasticFactory::registerMaterial<
    Marmot::Materials::LinearElasticInterface >( "LINEARELASTICINTERFACE" );

  (void)registered;
}
