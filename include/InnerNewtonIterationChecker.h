#pragma once
#include "bftTypedefs.h"

namespace bft{
    class InnerNewtonIterationChecker
    {
        const MatrixXd& residualScaleMatrix;
        const int nMaxInnerNewtonCycles;
        const int nMaxInnerNewtonCyclesAlt;
        const double innerNewtonTol;
        const double innerNewtonRTol;
        const double innerNewtonTolAlt;
        const double innerNewtonRTolAlt;

        public:
                
            InnerNewtonIterationChecker(const MatrixXd& residualScaleMatrix,
                                        int nMaxInnerNewtonCycles,
                                        int nMaxInnerNewtonCyclesAlt,
                                        double innerNewtonTol,
                                        double innerNewtonRTol,
                                        double innerNewtonTolAlt,
                                        double innerNewtonRTolAlt); 
            double relativeNorm(const VectorXd& increment, const VectorXd& reference);
            double residualNorm(const VectorXd& Residual);
            bool iterationFinished(const VectorXd& residual, const VectorXd& X, const VectorXd& dX, int numberOfIterations);
            bool isConverged(const VectorXd& residual, const VectorXd& X, const VectorXd& dX, int numberOfIterations);

    };
}
