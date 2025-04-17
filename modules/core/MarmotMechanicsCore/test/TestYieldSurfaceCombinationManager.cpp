#include "Marmot/MarmotTesting.h"
#include "Marmot/YieldSurfaceCombinationManager.h"

using namespace Marmot::Testing;

void testYieldSurfaceCombinationManager()
{
  Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< 3 >                   manager;
  Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< 3 >::YieldSurfFlagArr activeSurfaces;

  int                nCombinations = 0;
  std::ostringstream oss;

  while ( manager.getAnotherYieldFlagCombination( activeSurfaces ) ) {
    nCombinations++;
    manager.markYieldFlagCombinationAsUsed( activeSurfaces );

    oss << "Combination: ";
    for ( int i = 0; i < 3; i++ )
      oss << activeSurfaces( i ) << " ";
    oss << std::endl;
  }

  if ( nCombinations != 7 )
    throwExceptionOnFailure( false, "nCombinations != 7, but " + std::to_string( nCombinations ) );
}

int main()
{

  testYieldSurfaceCombinationManager();

  return 0;
}
