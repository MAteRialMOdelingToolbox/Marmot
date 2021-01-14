/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::NumericalAlgorithms {
    class InnerNewtonIterationChecker {

        const Eigen::MatrixXd& residualScaleMatrix;
        const int              nMaxInnerNewtonCycles;
        const int              nMaxInnerNewtonCyclesAlt;
        const double           innerNewtonTol;
        const double           innerNewtonRTol;
        const double           innerNewtonTolAlt;
        const double           innerNewtonRTolAlt;

      public:
        InnerNewtonIterationChecker( const Eigen::MatrixXd& residualScaleMatrix,
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
} // namespace Marmot::NumericalAlgorithms
