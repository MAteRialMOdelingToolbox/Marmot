#include "Marmot/InnerNewtonIterationChecker.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {

    InnerNewtonIterationChecker::InnerNewtonIterationChecker( const MatrixXd& residualScaleMatrix,
                                                              int             nMaxInnerNewtonCycles,
                                                              int             nMaxInnerNewtonCyclesAlt,
                                                              double          innerNewtonTol,
                                                              double          innerNewtonRTol,
                                                              double          innerNewtonTolAlt,
                                                              double          innerNewtonRTolAlt )
        : residualScaleMatrix( residualScaleMatrix ),
          nMaxInnerNewtonCycles( nMaxInnerNewtonCycles ),
          nMaxInnerNewtonCyclesAlt( nMaxInnerNewtonCyclesAlt ),
          innerNewtonTol( innerNewtonTol ),
          innerNewtonRTol( innerNewtonRTol ),
          innerNewtonTolAlt( innerNewtonTolAlt ),
          innerNewtonRTolAlt( innerNewtonRTolAlt )
    {
    }

    double InnerNewtonIterationChecker::relativeNorm( const VectorXd& increment, const VectorXd& reference )
    {
        double incNorm = increment.norm();
        double refNorm = reference.norm();

        if ( incNorm < 1e-16 )
            return incNorm;

        if ( refNorm < 1e-14 )
            // for a too small reference norm, a reasonable relative norm cannot be computed
            return 0.0;

        return increment.norm() / refNorm;
    }

    double InnerNewtonIterationChecker::residualNorm( const VectorXd& residual )
    {
        return ( residualScaleMatrix * residual ).norm();
    }

    bool InnerNewtonIterationChecker::iterationFinished( const VectorXd& residual,
                                                         const VectorXd& X,
                                                         const VectorXd& dX,
                                                         int             numberOfIterations )
    {
        if ( isConverged( residual, X, dX, numberOfIterations ) || numberOfIterations > nMaxInnerNewtonCyclesAlt )
            return true;
        else
            return false;
    }

    bool InnerNewtonIterationChecker::isConverged( const VectorXd& residual,
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
        else if ( numberOfIterations <= nMaxInnerNewtonCyclesAlt ) {
            if ( resNorm <= innerNewtonTolAlt && relNorm <= innerNewtonRTolAlt )
                return true;
        }
        return false;
    }
} // namespace Marmot
