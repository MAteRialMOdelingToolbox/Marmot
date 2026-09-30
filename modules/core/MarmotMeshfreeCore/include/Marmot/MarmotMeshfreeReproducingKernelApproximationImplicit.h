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
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::MarmotMeshfreeReproducingKernelApproximationImplicit
   * @brief Reproducing kernel approximation with implicit gradients.
   *
   * The shape function values are those of MarmotMeshfreeReproducingKernelApproximation. The gradients are not
   * obtained by differentiating @f$ \Psi_A @f$ but are constructed directly as corrected kernels (implicit, or
   * synchronized, gradients),
   * @f[
   *   \Psi^{(i)}_A(\boldsymbol{x}) = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A)\,
   *   \boldsymbol{b}^{(i)}(\boldsymbol{x})\, \phi_A(\boldsymbol{x}),
   *   \qquad \boldsymbol{M}(\boldsymbol{x})\, \boldsymbol{b}^{(i)}(\boldsymbol{x}) = \boldsymbol{H}^{(i)}_0,
   * @f]
   * with the same moment matrix @f$ \boldsymbol{M} @f$. The right-hand side @f$ \boldsymbol{H}^{(i)}_0 =
   * -\partial \boldsymbol{H}(\boldsymbol{z}) / \partial z_i |_{\boldsymbol{z} = \boldsymbol{0}} @f$ (the sign
   * stems from the shifted argument @f$ \boldsymbol{z} = \boldsymbol{x} - \boldsymbol{x}_A @f$) is @f$ -1 @f$ at
   * the position of the linear monomial @f$ z_i @f$ in @f$ \boldsymbol{H} @f$
   * and zero elsewhere (the position is found from Math::computeMonomialBasisGradient() at the origin, so it is
   * correct for any order). This enforces the gradient reproducing conditions
   * @f[
   *   \sum_A \Psi^{(i)}_A(\boldsymbol{x})\, \boldsymbol{x}_A^{\boldsymbol{\alpha}} =
   *   \frac{\partial \boldsymbol{x}^{\boldsymbol{\alpha}}}{\partial x_i}, \qquad |\boldsymbol{\alpha}| \le n,
   * @f]
   * so the implicit gradients reproduce the derivatives of all polynomials up to the completeness order, but they
   * are not the derivatives of the values @f$ \Psi_A @f$. They require neither kernel gradients nor
   * @f$ \boldsymbol{M}_{,i} @f$: all @f$ 1 + d @f$ right-hand sides are solved with one QR decomposition of
   * @f$ \boldsymbol{M} @f$.
   */
  class MarmotMeshfreeReproducingKernelApproximationImplicit : public MarmotMeshfreeReproducingKernelApproximation {

  public:
    /**
     * @brief Construct the approximation.
     * @param[in] dim               Spatial dimension @f$ d @f$.
     * @param[in] completenessOrder Desired completeness order @f$ n @f$.
     */
    MarmotMeshfreeReproducingKernelApproximationImplicit( int dim, int completenessOrder );

    /// @brief Default destructor.
    virtual ~MarmotMeshfreeReproducingKernelApproximationImplicit() = default;

    /**
     * @brief Not implemented; use computeShapeFunctionsAndGradients().
     * @param[in]  coord                       Evaluation point.
     * @param[in]  kernelFunctions             Candidate kernel functions.
     * @param[out] shapeFunctionValueGradients Not written.
     * @throws std::runtime_error always.
     */
    virtual void computeShapeFunctionGradients(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
      double*                                                   shapeFunctionValueGradients ) const override;

    /**
     * @brief Compute the RK shape function values @f$ \Psi_A @f$ and the implicit gradients @f$ \Psi^{(i)}_A @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     * @param[out] shapeFunctionValueGradients_ColMajor Implicit gradients, column-major
     *                                 @f$ d \times n_\mathrm{c} @f$ (entry @f$ i + d A @f$ is
     *                                 @f$ \Psi^{(i)}_A @f$); zero for non-covering candidates.
     * @throws std::runtime_error if the basis is empty (negative completeness order).
     */
    virtual void computeShapeFunctionsAndGradients(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
      double*                                                   shapeFunctionValues,
      double*                                                   shapeFunctionValueGradients_ColMajor ) const override;
  };

} // namespace Marmot::Meshfree
