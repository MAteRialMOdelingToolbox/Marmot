#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotTesting.h"
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Marmot::Testing;
using namespace MarmotLibrary;

// ---------------------------------------------------------------------------------------------
// createMaterial(): the success path is exercised elsewhere (e.g. CompressibleNeoHooke's own
// tests), but the throw branch for an unregistered material name is not. "COMPRESSIBLENEOHOOKE"
// is genuinely registered (by CompressibleNeoHookeRegistration.cpp, linked into every test
// executable via libMarmot), so a name guaranteed not to collide with it proves the throw path
// specifically, not just "any unknown name".
// ---------------------------------------------------------------------------------------------
void testCreateMaterialThrowsForUnregisteredName()
{
  bool threw = false;
  try {
    MarmotMaterialFiniteStrainFactory::createMaterial( "DEFINITELY_NOT_REGISTERED", nullptr, 0, 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "createMaterial() must throw std::invalid_argument for an unregistered material name "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testCreateMaterialSucceedsForRegisteredName()
{
  const double                                  properties[2] = { 3500., 1500. };
  std::unique_ptr< MarmotMaterialFiniteStrain > mat(
    MarmotMaterialFiniteStrainFactory::createMaterial( "COMPRESSIBLENEOHOOKE", properties, 2, 7 ) );

  throwExceptionOnFailure( mat != nullptr && mat->materialNumber == 7,
                           "createMaterial() did not return a correctly constructed instance in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testCreateMaterialThrowsForUnregisteredName,
    testCreateMaterialSucceedsForRegisteredName,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
