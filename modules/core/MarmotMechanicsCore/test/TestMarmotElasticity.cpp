#include "Marmot/MarmotTesting.h"

/*
 * ASSIGNEE: Manuel Hradsky
 * TODO: Test compliance and stiffness tensors for multiple cases and edge cases
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
