#pragma once
#include "bftTypedefs.h"

namespace bft{
    template <int nYieldSurfaces>
    class YieldSurfaceCombinationManager
    {
            const int idxUsedFlag;
            Array<bool, (2 << (nYieldSurfaces-1)), (nYieldSurfaces+1)> yieldSurfaceCombinations;

    public:
            YieldSurfaceCombinationManager();
            void initYieldFlagCombinations();
            bool getAnotherYieldFlagCombination(YieldSurfFlagArr<nYieldSurfaces>& activeSurfaces);
            void markYieldFlagCombinationAsUsed(const YieldSurfFlagArr<nYieldSurfaces>& activeSurfaces );
            void resetUsedYieldFlagCombinations();
    };
}

#include <cmath>

namespace bft{
    template <int n>
    YieldSurfaceCombinationManager<n>::YieldSurfaceCombinationManager():
            idxUsedFlag(n)
    {
        yieldSurfaceCombinations.setZero();
        const int numRows = 2<<(n-1);
        const int numCols = n+1;

        for(int row = 0; row < numRows; row++)
            for(int col = 0; col < numCols-1; col++)
                yieldSurfaceCombinations(row,col) = (row+1)& static_cast<int>(std::pow(2.0,col));

        yieldSurfaceCombinations.col(numCols-1).setConstant(false);
        return;
    }

    template <int n>
    bool YieldSurfaceCombinationManager<n>::getAnotherYieldFlagCombination(YieldSurfFlagArr<n>& activeSurfaces)
    {
        for(int i = 0; i < yieldSurfaceCombinations.rows(); i++)
        {

            bool alreadyUsed = yieldSurfaceCombinations.row(i)(idxUsedFlag);
            YieldSurfFlagArr<n> combinationCandidate = yieldSurfaceCombinations.row(i).head(n);
            if(!alreadyUsed && !(combinationCandidate == activeSurfaces).all() )
            {
                    activeSurfaces = combinationCandidate;
                    return true;
            }
        }
        return false;
    }

    template <int n>
    void YieldSurfaceCombinationManager<n>::resetUsedYieldFlagCombinations()
    {
            yieldSurfaceCombinations.col(idxUsedFlag).setConstant(false);
    }

    template <int n>
    void YieldSurfaceCombinationManager<n>::markYieldFlagCombinationAsUsed(const YieldSurfFlagArr<n>& activeSurfaces )
    {
            for(int i = 0; i < yieldSurfaceCombinations.rows(); i++)
                    if( (yieldSurfaceCombinations.row(i).head(n) == activeSurfaces).all())
                            yieldSurfaceCombinations.row(i)(idxUsedFlag) = true;
    }
}
