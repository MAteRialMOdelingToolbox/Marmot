#pragma once
#include "MarmotTypedefs.h"

/* Manager for yield surface combinations for multisurface plasticity
 * Matthias Neuner (2015)
 * */

namespace Marmot {
    template <int nYieldSurfaces>
    class YieldSurfaceCombinationManager {
        const int idxUsedFlag;

      public:
        Eigen::Array<bool, ( 1 << nYieldSurfaces ) - 1, ( nYieldSurfaces + 1 )> yieldSurfaceCombinations;
        typedef Eigen::Array<bool, 1, nYieldSurfaces>                           YieldSurfFlagArr;
        typedef Eigen::Array<double, 1, nYieldSurfaces>                         YieldSurfResArr;

        YieldSurfaceCombinationManager();
        void initYieldFlagCombinations();
        bool getAnotherYieldFlagCombination( YieldSurfFlagArr& activeSurfaces );
        void markYieldFlagCombinationAsUsed( const YieldSurfFlagArr& activeSurfaces );
        void resetUsedYieldFlagCombinations();
    };
} // namespace Marmot

namespace Marmot {
    template <int n>
    YieldSurfaceCombinationManager<n>::YieldSurfaceCombinationManager() : idxUsedFlag( n )
    {
        yieldSurfaceCombinations.setZero();
        const int numRows = ( 1 << n ) - 1;
        const int numCols = n;

        for ( int row = 0; row < numRows; row++ )
            for ( int col = 0; col < numCols; col++ )
                yieldSurfaceCombinations( row, col ) = ( row + 1 ) & 1 << ( col );
        return;
    }

    template <int n>
    bool YieldSurfaceCombinationManager<n>::getAnotherYieldFlagCombination( YieldSurfFlagArr& activeSurfaces )
    {
        for ( int i = 0; i < yieldSurfaceCombinations.rows(); i++ ) {
            bool             alreadyUsed          = yieldSurfaceCombinations.row( i )( idxUsedFlag );
            YieldSurfFlagArr combinationCandidate = yieldSurfaceCombinations.row( i ).head( n );
            if ( !alreadyUsed && ( combinationCandidate == true ).any() ) {
                activeSurfaces = combinationCandidate;
                return true;
            }
        }
        return false;
    }

    template <int n>
    void YieldSurfaceCombinationManager<n>::resetUsedYieldFlagCombinations()
    {
        yieldSurfaceCombinations.col( idxUsedFlag ).setConstant( false );
    }

    template <int n>
    void YieldSurfaceCombinationManager<n>::markYieldFlagCombinationAsUsed( const YieldSurfFlagArr& activeSurfaces )
    {
        for ( int i = 0; i < yieldSurfaceCombinations.rows(); i++ )
            if ( ( yieldSurfaceCombinations.row( i ).head( n ) == activeSurfaces ).all() )
                yieldSurfaceCombinations.row( i )( idxUsedFlag ) = true;
    }
} // namespace Marmot
