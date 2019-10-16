#pragma once
#include "bftTypedefs.h"

namespace bft {
    class InnerNewtonIterationCheckerMarkII{

        const Eigen::VectorXd residualScaleVector;
        const int              nMaxInnerNewtonCycles;
        const int              nMaxInnerNewtonCyclesAlt;
        const double           innerNewtonTol;
        const double           innerNewtonRTol;
        const double           innerNewtonTolAlt;
        const double           innerNewtonRTolAlt;

      public:
        InnerNewtonIterationCheckerMarkII( const Eigen::VectorXd& residualScaleVector,
                                     int                    nMaxInnerNewtonCycles,
                                     int                    nMaxInnerNewtonCyclesAlt,
                                     double                 innerNewtonTol,
                                     double                 innerNewtonRTol,
                                     double                 innerNewtonTolAlt,
                                     double                 innerNewtonRTolAlt );

        double relativeNorm( const Eigen::VectorXd& increment, const Eigen::VectorXd& reference );

        double residualNorm( const Eigen::VectorXd& Residual );

        bool iterationFinished( const Eigen::VectorXd& residual,
                                const Eigen::VectorXd& X,
                                const Eigen::VectorXd& dX,
                                int                    numberOfIterations );

        bool isConverged( const Eigen::VectorXd& residual,
                          const Eigen::VectorXd& X,
                          const Eigen::VectorXd& dX,
                          int                    numberOfIterations );
    };
} // namespace bft
