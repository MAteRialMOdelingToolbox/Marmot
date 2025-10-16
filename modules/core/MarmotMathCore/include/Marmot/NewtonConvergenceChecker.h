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
#include <Eigen/Core>

namespace Marmot::NumericalAlgorithms {

  /**
   * @class NewtonConvergenceChecker
   * @brief A class to check convergence of Newton-Raphson iterations.
   * This class provides methods to evaluate the convergence of Newton-Raphson iterations
   * based on residual norms and relative increments.
   */
  class NewtonConvergenceChecker {

    /// A vector to scale the residual for convergence checks.
    const Eigen::VectorXd residualScaleVector;
    /// Maximum number of Newton cycles allowed.
    const int nMaxNewtonCycles;
    /// Maximum number of alternative Newton cycles allowed.
    const int nMaxNewtonCyclesAlt;
    /// Tolerance for the Newton convergence based on absolute residual norm.
    const double newtonTol;
    /// Tolerance for the Newton convergence based on relative increment norm.
    const double newtonRTol;
    /// Alternative tolerance for the Newton convergence based on absolute residual norm.
    const double newtonTolAlt;
    /// Alternative tolerance for the Newton convergence based on relative increment norm.
    const double newtonRTolAlt;

  public:
    /**
     * @brief Constructor for NewtonConvergenceChecker.
     * @param residualScaleVector A vector to scale the residual for convergence checks.
     * @param nMaxNewtonCycles Maximum number of Newton cycles allowed.
     * @param nMaxNewtonCyclesAlt Maximum number of alternative Newton cycles allowed.
     * @param newtonTol Tolerance for the Newton convergence based on absolute residual norm.
     * @param newtonRTol Tolerance for the Newton convergence based on relative increment norm.
     * @param newtonTolAlt Alternative tolerance for the Newton convergence based on absolute residual norm.
     * @param newtonRTolAlt Alternative tolerance for the Newton convergence based on relative increment norm.
     */
    NewtonConvergenceChecker( const Eigen::VectorXd& residualScaleVector,
                              int                    nMaxNewtonCycles,
                              int                    nMaxNewtonCyclesAlt,
                              double                 newtonTol,
                              double                 newtonRTol,
                              double                 newtonTolAlt,
                              double                 newtonRTolAlt );

    /**
     * @brief Compute the relative norm of the increment with respect to a reference vector.
     * @param increment The increment vector.
     * @param reference The reference vector.
     *
     * The relative norm is computed as
     * \f[\displaystyle
     * r_{\rm rel} =
     * \begin{cases}
     * || \Delta \boldsymbol{X} ||_2 & \text{if } || \Delta \boldsymbol{X} ||_2 < 10^{-14} \\
     * 0.0 & \text{if } || \boldsymbol{X}_{\rm ref} ||_2 < 10^{-12} \\
     * \frac{ || \Delta \boldsymbol{X} ||_2 }{ || \boldsymbol{X}_{\rm ref} ||_2 } & \text{else}
     * \end{cases}
     * \f]
     * where \f$ || \square ||_2 \f$ denotes the Euclidean norm.
     *
     * @return The relative norm of the increment.
     */
    double relativeNorm( const Eigen::VectorXd& increment, const Eigen::VectorXd& reference );

    /**
     * @brief Compute the norm of the residual vector.
     * @param Residual The residual vector.
     * @details The residual is first scaled by the residualScaleVector before computing the norm.
     *
     * This function computes
     * \f[
     * r = \sqrt{ R^{\rm s}_i R^{\rm s}_i }  = || \boldsymbol{R}^{\rm s} ||_2
     * \f]
     * where \f$ R^{\rm s}_i = s_i R_i \f$ (no summation over i) are the scaled residual components
     * and \f$ s_i \f$ are the components of the residualScaleVector.
     *
     * @return The norm of the scaled residual.
     */
    double residualNorm( const Eigen::VectorXd& Residual );

    /**
     * @brief Check if the iteration has finished based on residuals and increments.
     * @param residual The current residual vector.
     * @param X The current solution vector.
     * @param dX The current increment vector.
     * @param numberOfIterations The current iteration count.
     * @return True if the iteration has finished, false otherwise.
     *
     * The iteration is considered finished if either the solution has converged
     * or the maximum number of alternative Newton cycles has been excee:ed:
     *```cpp
     * if ( isConverged( residual, X, dX, numberOfIterations ) || numberOfIterations > nMaxNewtonCyclesAlt )
     *   return true;
     * else
     *  return false;
     * ```
     */
    bool iterationFinished( const Eigen::VectorXd& residual,
                            const Eigen::VectorXd& X,
                            const Eigen::VectorXd& dX,
                            int                    numberOfIterations );
    /**
     * @brief Check if the solution has converged based on residuals and increments.
     * @param residual The current residual vector.
     * @param X The current solution vector.
     * @param dX The current increment vector.
     * @param numberOfIterations The current iteration count.
     * @return True if the solution has converged, false otherwise.
     *
     * The solution is considered converged if both the absolute residual norm
     * and the relative increment norm is below the specified tolerances.
     * Alternative tolerances are used if the number of iterations exceeds nMaxNewtonCycles.
     *
     * ```cpp
     * if ( numberOfIterations <= nMaxNewtonCycles ) {
     *   if ( relativeNorm( dX, X ) <= newtonRTol && residualNorm( residual ) <= newtonTol )
     *     return true;
     * } else if ( numberOfIterations <= nMaxNewtonCyclesAlt + 1 ){
     *   if ( relativeNorm( dX, X ) <= newtonRTolAlt && residualNorm( residual ) <= newtonTolAlt )
     *     return true;
     * }
     * else
     *   return false;
     * ```
     */
    bool isConverged( const Eigen::VectorXd& residual,
                      const Eigen::VectorXd& X,
                      const Eigen::VectorXd& dX,
                      int                    numberOfIterations );
  };
} // namespace Marmot::NumericalAlgorithms
