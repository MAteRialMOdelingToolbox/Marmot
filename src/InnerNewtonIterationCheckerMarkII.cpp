#include "InnerNewtonIterationCheckerMarkII.h"
#include <iostream>

using namespace Eigen;

namespace bft {

    InnerNewtonIterationCheckerMarkII::InnerNewtonIterationCheckerMarkII( const VectorXd& residualScaleVector,
                                                              int             nMaxInnerNewtonCycles,
                                                              int             nMaxInnerNewtonCyclesAlt,
                                                              double          innerNewtonTol,
                                                              double          innerNewtonRTol,
                                                              double          innerNewtonTolAlt,
                                                              double          innerNewtonRTolAlt )
        : residualScaleVector( residualScaleVector ),
          nMaxInnerNewtonCycles( nMaxInnerNewtonCycles ),
          nMaxInnerNewtonCyclesAlt( nMaxInnerNewtonCyclesAlt ),
          innerNewtonTol( innerNewtonTol ),
          innerNewtonRTol( innerNewtonRTol ),
          innerNewtonTolAlt( innerNewtonTolAlt ),
          innerNewtonRTolAlt( innerNewtonRTolAlt )
    {
    }

    double InnerNewtonIterationCheckerMarkII::relativeNorm( const VectorXd& increment, const VectorXd& reference )
    {
        double incNorm = increment.norm();
        double refNorm = reference.norm();

        if ( incNorm < 1e-14 )
            return incNorm;

        if ( refNorm < 1e-12 )
            // for a too small reference norm, a reasonable relative norm cannot be computed
            return 0.0;

        return increment.norm() / refNorm;
    }

    double InnerNewtonIterationCheckerMarkII::residualNorm( const VectorXd& residual )
    {
        return ( residualScaleVector.array()* residual.array() ).matrix().norm();
    }

    bool InnerNewtonIterationCheckerMarkII::iterationFinished( const VectorXd& residual,
                                                         const VectorXd& X,
                                                         const VectorXd& dX,
                                                         int             numberOfIterations )
    {
        if ( isConverged( residual, X, dX, numberOfIterations ) || numberOfIterations > nMaxInnerNewtonCyclesAlt )
            return true;
        else
            return false;
    }

    bool InnerNewtonIterationCheckerMarkII::isConverged( const VectorXd& residual,
                                                   const VectorXd& X,
                                                   const VectorXd& dX,
                                                   int             numberOfIterations )
    {
        const double resNorm = residualNorm( residual );
        const double relNorm = relativeNorm( dX, X );

        if ( numberOfIterations <= nMaxInnerNewtonCycles ) {
            if ( resNorm <= innerNewtonTol && relNorm <= innerNewtonRTol )
                return true;
        }
        else if ( numberOfIterations <= nMaxInnerNewtonCyclesAlt +1 ) {
            if ( resNorm <= innerNewtonTolAlt && relNorm <= innerNewtonRTolAlt )
                return true;
        }
        return false;
    }
} // namespace bft
