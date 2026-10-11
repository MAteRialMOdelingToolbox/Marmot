#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement::Quadrature::Spatial2D;

// ---------------------------------------------------------------------------------------------
// modifyCharElemLengthAbaqusLike(): has no callers anywhere in the codebase (verified by search),
// so it is tested directly against its own documented per-integration-point scale factors for a
// Quad8-style 3x3 (9-point) integration rule.
// ---------------------------------------------------------------------------------------------

void testModifyCharElemLengthAbaqusLikeCentralNode()
{
  double length = 1.0;
  modifyCharElemLengthAbaqusLike( length, 0 );
  throwExceptionOnFailure( checkIfEqual( length, 2.0 / 3.0, 1e-12 ),
                           "The central integration point (0) must scale the length by 2/3." );
}

void testModifyCharElemLengthAbaqusLikeCornerNodes()
{
  for ( int intPoint = 1; intPoint <= 4; intPoint++ ) {
    double length = 1.0;
    modifyCharElemLengthAbaqusLike( length, intPoint );
    throwExceptionOnFailure( checkIfEqual( length, 5.0 / 12.0, 1e-12 ),
                             "Corner integration points (1-4) must scale the length by 5/12." );
  }
}

void testModifyCharElemLengthAbaqusLikeMiddleNodes()
{
  for ( int intPoint = 5; intPoint <= 8; intPoint++ ) {
    double length = 1.0;
    modifyCharElemLengthAbaqusLike( length, intPoint );
    throwExceptionOnFailure( checkIfEqual( length, std::sqrt( 5.0 / 18.0 ), 1e-12 ),
                             "Middle integration points (5-8) must scale the length by sqrt(5/18)." );
  }
}

void testModifyCharElemLengthAbaqusLikeLeavesLengthUnchangedForUnhandledIndex()
{
  // No default case: an out-of-range intPoint index is a documented-by-absence no-op.
  double length = 1.0;
  modifyCharElemLengthAbaqusLike( length, 99 );
  throwExceptionOnFailure( checkIfEqual( length, 1.0, 1e-12 ),
                           "An unhandled integration point index must leave the length unchanged." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testModifyCharElemLengthAbaqusLikeCentralNode,
    testModifyCharElemLengthAbaqusLikeCornerNodes,
    testModifyCharElemLengthAbaqusLikeMiddleNodes,
    testModifyCharElemLengthAbaqusLikeLeavesLengthUnchangedForUnhandledIndex,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
