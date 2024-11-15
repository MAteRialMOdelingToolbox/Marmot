#include "Marmot/MarmotTesting.h"

/*
 * ASSIGNEE: Nasser Alkmim
 * TODO: Test MarMotStateVarVectorManager implementation
 */

using namespace Marmot::Testing;

void myDummyTestFunction()
{
  throwExceptionOnFailure( checkIfEqual( 1, 1 ), "1 == 1" );
}

int main()
{
  myDummyTestFunction();
  return 0;
}
