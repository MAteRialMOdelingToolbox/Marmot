#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

class DummyMaterialStateVarManager : public MarmotStateVarVectorManager {
  /*
   * DummyMaterialStateVarManager is a derived class of MarmotStateVarVectorManager
   * for testing purposes. It defines a static layout with two state variables:
   * "var1" and "var2", with lengths 3 and 2 respectively.
   */

public:
  inline const static auto layout = makeLayout( { { .name = "var1", .length = 3 }, { .name = "var2", .length = 2 } } );

  // Constructor initializes the parent class
  DummyMaterialStateVarManager( double* stateVars ) : MarmotStateVarVectorManager( stateVars, layout ){};

  // Expose nRequiredStateVars for testing layout calculation
  static int getRequiredSize() { return layout.nRequiredStateVars; }
};

void testStateVarVectorManagerFind()
{
  /*
   * Test find method in the statevar vector manager
   * and check if the pointers to the state variables are correct.
   */
  std::vector< double >        stateVarData( 5 );
  DummyMaterialStateVarManager manager( stateVarData.data() );

  double& var1 = manager.find( "var1" );
  throwExceptionOnFailure( checkIfEqual( var1, 0.0 ), "Initial value of var1" );
  throwExceptionOnFailure( &var1 == stateVarData.data(), "var1 should point to the first element of stateVarData" );
  double& var2 = manager.find( "var2" );
  throwExceptionOnFailure( &var2 == stateVarData.data() + 3,
                           "var2 should point to the fourth element of stateVarData" );
}

void testStateVarVectorManagerContains()
{
  /*
   * Test contains method in the statevar vector manager
   * and check if the correct state variables are contained.
   */
  std::vector< double >        stateVarData( 5 );
  DummyMaterialStateVarManager manager( stateVarData.data() );

  throwExceptionOnFailure( manager.contains( "var1" ), "var1 should be contained" );
  throwExceptionOnFailure( manager.contains( "var2" ), "var2 should be contained" );
  throwExceptionOnFailure( !manager.contains( "var3" ), "var3 should not be contained" );
  throwExceptionOnFailure( !manager.contains( "" ), "Empty string should not be contained" );
}

void testStateVarVectorManagerLayout()
{
  /*
   * Test the layout of the statevar vector manager
   * and check if the required size is correct.
   */
  throwExceptionOnFailure( checkIfEqual( DummyMaterialStateVarManager::getRequiredSize(), 5 ),
                           "Layout requires 5 doubles" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{
    testStateVarVectorManagerFind,
    testStateVarVectorManagerContains,
    testStateVarVectorManagerLayout,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
