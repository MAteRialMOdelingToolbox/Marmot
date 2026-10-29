#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace MarmotLibrary;

namespace {

  // Minimal, registerable material used only to exercise the factory's registration/creation
  // machinery for nNonlocalVariables values (2 and 6) that no shipped material currently uses --
  // the factory itself is a general template over nNonlocalVariables, explicitly instantiated for
  // 1, 2 and 6, so this is genuine, reachable library code, not dead code.
  template < int nNonlocalVariables >
  class DummyMaterial : public MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables > {
  public:
    using Base = MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >;
    using Base::Base;

    void computeStress( typename Base::response&,
                        typename Base::tangents&,
                        const typename Base::increment& ) const override
    {
    }

    double getDensity( const double* ) const override { return 1.0; }

    std::vector< double > getNonlocalViscosity( const double* ) const override
    {
      return std::vector< double >( nNonlocalVariables, 0.0 );
    }
  };

} // namespace

// ---------------------------------------------------------------------------------------------
// createMaterial<1>(): the throw branch (unregistered material name) is never exercised by
// production code, since every caller only ever requests materials it just confirmed are
// registered. "AT2PHASEFIELD" is genuinely registered (by AT2PhaseFieldRegistration.cpp, linked
// into every test executable via libMarmot), so a name guaranteed not to collide with it proves
// the throw path specifically, not just "any unknown name".
// ---------------------------------------------------------------------------------------------
void testCreateMaterialN1ThrowsForUnregisteredName()
{
  bool threw = false;
  try {
    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 1 >::createMaterial( "DEFINITELY_NOT_REGISTERED",
                                                                                  nullptr,
                                                                                  0,
                                                                                  1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "createMaterial<1>() must throw std::invalid_argument for an unregistered material name "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// registerMaterial<2>() / createMaterial<2>(): round-trip registration and creation for the
// nNonlocalVariables=2 specialization.
// ---------------------------------------------------------------------------------------------
void testRegisterAndCreateMaterialN2()
{
  MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 2 >::registerMaterial< DummyMaterial< 2 > >(
    "TESTDUMMYMATERIALN2" );

  std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< 2 > > mat(
    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 2 >::createMaterial( "TESTDUMMYMATERIALN2",
                                                                                  nullptr,
                                                                                  0,
                                                                                  7 ) );

  throwExceptionOnFailure( mat != nullptr && mat->materialNumber == 7,
                           "createMaterial<2>() did not return a correctly constructed instance in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testCreateMaterialN2ThrowsForUnregisteredName()
{
  bool threw = false;
  try {
    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 2 >::createMaterial( "SOME_OTHER_UNKNOWN_NAME",
                                                                                  nullptr,
                                                                                  0,
                                                                                  1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "createMaterial<2>() must throw std::invalid_argument for an unregistered material name "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// registerMaterial<6>() / createMaterial<6>(): round-trip registration and creation for the
// nNonlocalVariables=6 specialization.
// ---------------------------------------------------------------------------------------------
void testRegisterAndCreateMaterialN6()
{
  MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::registerMaterial< DummyMaterial< 6 > >(
    "TESTDUMMYMATERIALN6" );

  std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 > > mat(
    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::createMaterial( "TESTDUMMYMATERIALN6",
                                                                                  nullptr,
                                                                                  0,
                                                                                  3 ) );

  throwExceptionOnFailure( mat != nullptr && mat->materialNumber == 3,
                           "createMaterial<6>() did not return a correctly constructed instance in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testCreateMaterialN6ThrowsForUnregisteredName()
{
  bool threw = false;
  try {
    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< 6 >::createMaterial( "SOME_OTHER_UNKNOWN_NAME",
                                                                                  nullptr,
                                                                                  0,
                                                                                  1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "createMaterial<6>() must throw std::invalid_argument for an unregistered material name "
                           "in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
int main()
{
  const std::vector< std::function< void() > > tests = {
    testCreateMaterialN1ThrowsForUnregisteredName,
    testRegisterAndCreateMaterialN2,
    testCreateMaterialN2ThrowsForUnregisteredName,
    testRegisterAndCreateMaterialN6,
    testCreateMaterialN6ThrowsForUnregisteredName,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
