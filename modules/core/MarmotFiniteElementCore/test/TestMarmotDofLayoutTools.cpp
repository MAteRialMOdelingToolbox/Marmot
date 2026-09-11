#include "Marmot/MarmotDofLayoutTools.h"
#include "Marmot/MarmotTesting.h"
#include <algorithm>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement;

// ---------------------------------------------------------------------------------------------
// These utilities have no callers anywhere in the codebase (verified by search), so there is no
// existing usage to infer the intended semantics from. The pair<int,int> is (dofsPerNode,
// nodesForField): this matches how the implementation actually destructures and uses it -- the
// header's own doc comment had it backwards, and has been corrected alongside these tests.
// ---------------------------------------------------------------------------------------------

void testMakeNodeFieldLayoutSingleFieldOnAllNodes()
{
  const std::map< std::string, std::pair< int, int > > fieldSizes = { { "displacement", { 2, 4 } } };

  const auto nodeFields = makeNodeFieldLayout( fieldSizes );

  throwExceptionOnFailure( nodeFields.size() == 4, "makeNodeFieldLayout() returned the wrong number of nodes." );
  for ( const auto& fields : nodeFields ) {
    throwExceptionOnFailure( fields.size() == 1, "Every node should carry exactly one field." );
    throwExceptionOnFailure( fields[0] == "displacement", "The single field must be \"displacement\"." );
  }
}

void testMakeNodeFieldLayoutMixedOrderElement()
{
  // A Taylor-Hood-style mixed element: quadratic displacement (2 dofs/node) on 8 nodes, linear
  // pressure (1 dof/node) on only the first 4 (corner) nodes.
  const std::map< std::string, std::pair< int, int > > fieldSizes = { { "displacement", { 2, 8 } },
                                                                      { "pressure", { 1, 4 } } };

  const auto nodeFields = makeNodeFieldLayout( fieldSizes );

  throwExceptionOnFailure( nodeFields.size() == 8,
                           "makeNodeFieldLayout() must use the larger of the two fields' node counts." );

  for ( int i = 0; i < 4; i++ ) {
    throwExceptionOnFailure( nodeFields[i].size() == 2, "Corner nodes must carry both displacement and pressure." );
    throwExceptionOnFailure( nodeFields[i][0] == "displacement" && nodeFields[i][1] == "pressure",
                             "Field order must follow the (alphabetical) map iteration order." );
  }
  for ( int i = 4; i < 8; i++ ) {
    throwExceptionOnFailure( nodeFields[i].size() == 1 && nodeFields[i][0] == "displacement",
                             "Midside nodes must carry only displacement." );
  }
}

void testMakeBlockedLayoutPermutationPatternGroupsDofsByField()
{
  // Same mixed element as above: 8 nodes x 2 displacement dofs (16) + 4 nodes x 1 pressure dof
  // (4) = 20 total DOFs, laid out node-major as
  // [u0x,u0y,p0, u1x,u1y,p1, u2x,u2y,p2, u3x,u3y,p3, u4x,u4y, u5x,u5y, u6x,u6y, u7x,u7y].
  const std::map< std::string, std::pair< int, int > > fieldSizes = { { "displacement", { 2, 8 } },
                                                                      { "pressure", { 1, 4 } } };
  const auto                                           nodeFields = makeNodeFieldLayout( fieldSizes );

  const auto pattern = makeBlockedLayoutPermutationPattern( nodeFields, fieldSizes );

  const std::vector< int > expected = { 0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 2, 5, 8, 11 };

  throwExceptionOnFailure( pattern.size() == expected.size(),
                           "makeBlockedLayoutPermutationPattern() returned the wrong number of entries." );
  for ( size_t i = 0; i < expected.size(); i++ )
    throwExceptionOnFailure( pattern[i] == expected[i],
                             "makeBlockedLayoutPermutationPattern()[" + std::to_string( i ) + "] mismatch." );

  // The pattern must be a permutation of 0..19 (every DOF appears exactly once).
  std::vector< int > sorted = pattern;
  std::sort( sorted.begin(), sorted.end() );
  for ( int i = 0; i < 20; i++ )
    throwExceptionOnFailure( sorted[i] == i, "The permutation pattern must cover every DOF exactly once." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testMakeNodeFieldLayoutSingleFieldOnAllNodes,
    testMakeNodeFieldLayoutMixedOrderElement,
    testMakeBlockedLayoutPermutationPatternGroupsDofsByField,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
