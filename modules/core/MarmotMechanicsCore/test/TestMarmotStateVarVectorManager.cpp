#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

// Create dummy manager for testing (based on VonMises.h)
class DummyMaterialStateVarManager : public MarmotStateVarVectorManager {
public:
  inline const static auto layout = makeLayout( { { .name = "var1", .length = 3 }, { .name = "var2", .length = 2 } } );

  // Constructor initializes the parent class
  DummyMaterialStateVarManager( double* stateVars ) : MarmotStateVarVectorManager( stateVars, layout ) {};
};

void testStateVarVectorManagerFind()
{
  DummyMaterialStateVarManager manager( new double[5] );

  double& var1 = manager.find( "var1" );
  var1         = 1.0;
  throwExceptionOnFailure( checkIfEqual( 1, var1 ), "1 == 1" );
}

void testStateVarVectorManagerContains()
{
  DummyMaterialStateVarManager manager( new double[5] );

  throwExceptionOnFailure( manager.contains( "var1" ), "var1 should be contained" );
  throwExceptionOnFailure( manager.contains( "var2" ), "var2 should be contained" );
  throwExceptionOnFailure( !manager.contains( "var3" ), "var3 should not be contained" );
}

int main()
{
  testStateVarVectorManagerFind();
  testStateVarVectorManagerContains();

  return 0;
}
