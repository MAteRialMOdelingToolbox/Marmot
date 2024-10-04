#include <Marmot/MarmotDofLayoutTools.h>
#include <vector>

namespace Marmot {

  namespace FiniteElement {

    using namespace std;

    const vector< vector< string > > makeNodeFieldLayout( const map< string, pair< int, int > >& fieldSizes )

    {
      vector< vector< string > > nodeFields;

      int maxNumberOfNodes = 0;
      for ( const auto& [field, fieldSize] : fieldSizes ) {
        const auto& [nFieldComponents, nNodesForField] = fieldSize;
        maxNumberOfNodes                               = max( maxNumberOfNodes, nNodesForField );
      }

      for ( int idxNode = 0; idxNode < maxNumberOfNodes; idxNode++ ) {
        nodeFields.push_back( vector< string >() );

        for ( const auto& [field, fieldSize] : fieldSizes ) {
          const auto& [nFieldComponents, nNodesForField] = fieldSize;
          if ( idxNode < nNodesForField )
            nodeFields[idxNode].push_back( field );
        }
      }

      return nodeFields;
    }

    vector< int > makeBlockedLayoutPermutationPattern( const vector< vector< string > >&      nodeFields,
                                                       const map< string, pair< int, int > >& fieldSizes )
    {

      vector< int > permutationPattern;

      for ( const auto& [blockedField, fieldSize] : fieldSizes ) {

        int currentIdx = 0;
        for ( size_t i = 0; i < nodeFields.size(); i++ ) {
          for ( size_t j = 0; j < nodeFields[i].size(); j++ ) {

            const auto& currentNodeField                         = nodeFields[i][j];
            const auto [nCurrentFieldComponents, nNodesForField] = fieldSizes.at( currentNodeField );

            if ( blockedField == currentNodeField )
              for ( int component = 0; component < nCurrentFieldComponents; component++ )
                permutationPattern.push_back( currentIdx + component );

            currentIdx += nCurrentFieldComponents;
          }
        }
      }

      return permutationPattern;
    }

  } // namespace FiniteElement
} // namespace Marmot
