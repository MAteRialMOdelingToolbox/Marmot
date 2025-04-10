#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;

// Create dummy manager for testing (based on VonMises.h)
// Base class is protected, so we need to create a derived class
class DummyMaterialStateVarManager : public MarmotStateVarVectorManager {
public:
  inline const static auto layout = makeLayout( {
    { .name = "var1", .length = 3 },
    { .name = "var2", .length = 2 } } );

  // Constructor initializes the parent class
  DummyMaterialStateVarManager( double* stateVars ) : MarmotStateVarVectorManager( stateVars, layout ) {};
};

void testStateVarVectorManagerFind()
{
  std::vector<double> stateVarData(5);
  DummyMaterialStateVarManager manager( stateVarData.data() );

  double& var1 = manager.find( "var1" );
  throwExceptionOnFailure( checkIfEqual( var1, 0.0 ), "Initial value of var1" );
  throwExceptionOnFailure(&var1 == stateVarData.data(), "var1 should point to the first element of stateVarData" );
  double& var2 = manager.find( "var2" );
  throwExceptionOnFailure(&var2 == stateVarData.data() + 3, "var2 should point to the fourth element of stateVarData" );
}

void testStateVarVectorManagerContains()
{
  std::vector<double> stateVarData(5);
  DummyMaterialStateVarManager manager( stateVarData.data() );

  throwExceptionOnFailure( manager.contains( "var1" ), "var1 should be contained" );
  throwExceptionOnFailure( manager.contains( "var2" ), "var2 should be contained" );
  throwExceptionOnFailure( !manager.contains( "var3" ), "var3 should not be contained" );
  throwExceptionOnFailure( !manager.contains( "" ), "Empty string should not be contained" );
}

int main()
{
  testStateVarVectorManagerFind();
  testStateVarVectorManagerContains();

  return 0;
}
