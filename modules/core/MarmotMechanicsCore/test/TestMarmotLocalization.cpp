#include "Marmot/MarmotTesting.h"

/*
 * ASSIGNEE: Paul Hofer
 * TODO: Test MarmotLocalization for multiple cases and edge cases
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
