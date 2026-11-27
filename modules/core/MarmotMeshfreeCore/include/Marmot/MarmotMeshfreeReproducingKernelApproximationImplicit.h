/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
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
#include <limits>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::ReproducingKernelApproximationImplicit
   * @brief Reproducing kernel approximation with implicit gradients.
   *
   * The shape function values are those of ReproducingKernelApproximation. The gradients are not
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
   * and zero elsewhere (the position is found from the gradient of the basis at the origin, Math::MonomialBasis, so it
   * is correct for any order). This enforces the gradient reproducing conditions
   * @f[
   *   \sum_A \Psi^{(i)}_A(\boldsymbol{x})\, \boldsymbol{x}_A^{\boldsymbol{\alpha}} =
   *   \frac{\partial \boldsymbol{x}^{\boldsymbol{\alpha}}}{\partial x_i}, \qquad |\boldsymbol{\alpha}| \le n,
   * @f]
   * so the implicit gradients reproduce the derivatives of all polynomials up to the completeness order, but they
   * are not the derivatives of the values @f$ \Psi_A @f$. They require neither kernel gradients nor
   * @f$ \boldsymbol{M}_{,i} @f$: all @f$ 1 + d @f$ right-hand sides are solved with one QR decomposition of
   * @f$ \boldsymbol{M} @f$. The completeness order of a point is chosen as for ReproducingKernelApproximation,
   * including the optional limit on the condition of @f$ \boldsymbol{M} @f$.
   *
   * At a point evaluated with the order 0 (too few covering kernels, or the limit on the condition), the implicit
   * gradients would reproduce only the derivative of a constant and vanish identically: such a point would carry
   * no strain. It gets the Shepard gradients instead, the derivatives of its values
   * @f$ \Psi_A = \phi_A / \sum_B \phi_B @f$, as with explicit gradients; only there are kernel gradients evaluated.
   */
  class ReproducingKernelApproximationImplicit : public ReproducingKernelApproximation {

  public:
    /**
     * @brief Construct the approximation.
     * @param[in] dim                        Spatial dimension @f$ d \le @f$ maxDim.
     * @param[in] completenessOrder          Desired completeness order @f$ n \ge 0 @f$.
     * @param[in] momentMatrixConditionLimit The limit @f$ \kappa_\mathrm{max} \ge 1 @f$ on the condition of the
     *                                       scaled moment matrix, above which the completeness order of a point
     *                                       is lowered; infinite (the default) for no limit (see
     *                                       ReproducingKernelApproximation).
     * @throws std::invalid_argument if the order is negative, the dimension or the basis exceeds maxDim or
     * maxSizeOfBasis, or the limit is below 1 or not a number (see ReproducingKernelApproximation).
     */
    ReproducingKernelApproximationImplicit(
      int    dim,
      int    completenessOrder,
      double momentMatrixConditionLimit = std::numeric_limits< double >::infinity() );

    /// @brief Default destructor.
    virtual ~ReproducingKernelApproximationImplicit() = default;

    /**
     * @brief Compute the RK shape function values @f$ \Psi_A @f$ and the implicit gradients @f$ \Psi^{(i)}_A @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     * @param[out] shapeFunctionValueGradients_ColMajor Implicit gradients, column-major
     *                                 @f$ d \times n_\mathrm{c} @f$ (entry @f$ i + d A @f$ is
     *                                 @f$ \Psi^{(i)}_A @f$; the Shepard gradients at a point of order 0); zero for
     *                                 non-covering candidates.
     * @throws std::runtime_error if the moment matrix is singular (see factorizeMomentMatrix()).
     */
    virtual void computeShapeFunctionsAndGradients( const double*                               coord,
                                                    const std::vector< const KernelFunction* >& kernelFunctions,
                                                    double*                                     shapeFunctionValues,
                                                    double* shapeFunctionValueGradients_ColMajor ) const override;
  };

} // namespace Marmot::Discretizations::Meshfree
