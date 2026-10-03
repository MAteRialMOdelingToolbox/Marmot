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

#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::MarmotMeshfreeApproximation
   * @brief Abstract interface of a meshfree approximation, which evaluates shape functions from a set of kernels.
   *
   * Given an evaluation point @f$ \boldsymbol{x} @f$ and a list of @f$ n @f$ candidate kernel functions
   * @f$ \phi_A @f$, @f$ A = 0, \dots, n-1 @f$ (MarmotMeshfreeKernelFunction), an approximation computes the
   * shape functions @f$ \Psi_A(\boldsymbol{x}) @f$ and/or their gradients. The output arrays are indexed by the
   * position of the kernel in the candidate list; candidates that do not cover @f$ \boldsymbol{x} @f$ get zero
   * entries.
   *
   * The gradient output is a column-major @f$ d \times n @f$ array, i.e., entry @f$ i + d\,A @f$ holds
   * @f$ \partial \Psi_A / \partial x_i @f$, where @f$ d @f$ is the spatial dimension.
   */
  class MarmotMeshfreeApproximation {

  public:
    /// @brief Default constructor.
    MarmotMeshfreeApproximation() = default;

    /// @brief Virtual destructor; approximations are used through base-class pointers.
    virtual ~MarmotMeshfreeApproximation() = default;

    /**
     * @brief Compute the shape function values @f$ \Psi_A(\boldsymbol{x}) @f$.
     * @param[in]  coord               Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n @f$).
     */
    virtual void computeShapeFunctions( const double*                                             coord,
                                        const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
                                        double* shapeFunctionValues ) const = 0;

    /**
     * @brief Compute the shape function gradients @f$ \partial \Psi_A / \partial x_i @f$.
     * @param[in]  coord                       Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point.
     * @param[in]  kernelFunctions             The @f$ n @f$ candidate kernel functions.
     * @param[out] shapeFunctionValueGradients Shape function gradients, column-major @f$ d \times n @f$.
     */
    virtual void computeShapeFunctionGradients(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
      double*                                                   shapeFunctionValueGradients ) const = 0;

    /**
     * @brief Compute the shape function values and gradients in one pass.
     * @param[in]  coord                       Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point.
     * @param[in]  kernelFunctions             The @f$ n @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues         Shape function values (length @f$ n @f$).
     * @param[out] shapeFunctionValueGradients Shape function gradients, column-major @f$ d \times n @f$.
     */
    virtual void computeShapeFunctionsAndGradients(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
      double*                                                   shapeFunctionValues,
      double*                                                   shapeFunctionValueGradients ) const = 0;
  };

} // namespace Marmot::Meshfree
