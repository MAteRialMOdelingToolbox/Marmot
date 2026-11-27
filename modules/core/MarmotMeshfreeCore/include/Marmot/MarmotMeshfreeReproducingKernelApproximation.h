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

#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/QR>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @class Marmot::Discretizations::Meshfree::ReproducingKernelApproximation
   * @brief Reproducing kernel (RK) approximation with explicit (direct) shape function gradients.
   *
   * The shape function of node @f$ A @f$ is the kernel @f$ \phi_A @f$ (KernelFunction) multiplied
   * by a correction function,
   * @f[
   *   \Psi_A(\boldsymbol{x}) = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A)\, \boldsymbol{b}(\boldsymbol{x})\,
   *   \phi_A(\boldsymbol{x}),
   * @f]
   * where @f$ \boldsymbol{H} @f$ is the vector of all monomials up to the completeness order @f$ n @f$
   * (Math::MonomialBasis, first entry @f$ 1 @f$) evaluated at the shifted coordinates
   * @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$. The coefficient vector @f$ \boldsymbol{b} @f$ solves
   * @f[
   *   \boldsymbol{M}(\boldsymbol{x})\, \boldsymbol{b}(\boldsymbol{x}) = \boldsymbol{H}_0,
   *   \qquad \boldsymbol{M}(\boldsymbol{x}) = \sum_{B} \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_B)\,
   *   \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_B)\, \phi_B(\boldsymbol{x}),
   *   \qquad \boldsymbol{H}_0 = \boldsymbol{H}(\boldsymbol{0}) = [1, 0, \dots, 0]^T,
   * @f]
   * with the moment matrix @f$ \boldsymbol{M} @f$ summed over the kernels covering @f$ \boldsymbol{x} @f$. This
   * enforces the reproducing conditions
   * @f$ \sum_A \Psi_A(\boldsymbol{x})\, (\boldsymbol{x} - \boldsymbol{x}_A)^{\boldsymbol{\alpha}} =
   * \delta_{\boldsymbol{\alpha}\boldsymbol{0}} @f$, or equivalently
   * @f$ \sum_A \Psi_A(\boldsymbol{x})\, \boldsymbol{x}_A^{\boldsymbol{\alpha}} = \boldsymbol{x}^{\boldsymbol{\alpha}}
   * @f$, for all multi-indices @f$ |\boldsymbol{\alpha}| \le n @f$; in particular, the shape functions form a
   * partition of unity. The continuity of @f$ \Psi_A @f$ is that of the kernels.
   *
   * The gradients are the exact derivatives of @f$ \Psi_A @f$ (product rule),
   * @f[
   *   \Psi_{A,i} = \left( \boldsymbol{b}_{,i}^T \boldsymbol{H} + \boldsymbol{b}^T \boldsymbol{H}_{,i} \right)
   *   \phi_A + \boldsymbol{b}^T \boldsymbol{H}\, \phi_{A,i},
   *   \qquad \boldsymbol{b}_{,i} = -\boldsymbol{M}^{-1} \boldsymbol{M}_{,i}\, \boldsymbol{b},
   *   \qquad \boldsymbol{M}_{,i} = \sum_B \left( \boldsymbol{H}_{,i} \boldsymbol{H}^T + \boldsymbol{H}
   *   \boldsymbol{H}_{,i}^T \right) \phi_B + \boldsymbol{H} \boldsymbol{H}^T \phi_{B,i},
   * @f]
   * with @f$ (\bullet)_{,i} = \partial (\bullet) / \partial x_i @f$ and @f$ \boldsymbol{H} @f$ evaluated at
   * @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$ (resp. @f$ \boldsymbol{x} - \boldsymbol{x}_B @f$). The linear
   * systems are solved with a column-pivoting Householder QR decomposition of @f$ \boldsymbol{M} @f$, which is
   * computed once per evaluation point.
   *
   * **One evaluation, four passes over the kernels.** First, every candidate kernel is evaluated once, with its
   * gradient if the gradients of the shape functions are asked for (evaluateKernelFunctions()); the values and
   * gradients are kept in the output arrays, which have one entry per candidate. Second, the basis
   * @f$ \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_A) @f$ of every covering kernel (and its gradient) is evaluated
   * once and kept (evaluateBasisOfKernels()): the moment matrix and the shape functions both need it. Third,
   * @f$ \boldsymbol{M} @f$ (and @f$ \boldsymbol{M}_{,i} @f$) is assembled (assembleMomentMatrix()); fourth, the shape
   * functions follow, overwriting the kernel values. The moment matrix and the vectors of the basis are of a fixed
   * largest size (BasisVector, MomentMatrix); the bases of the kernels, one column per candidate, are the one
   * allocation of an evaluation.
   *
   * If fewer kernels cover @f$ \boldsymbol{x} @f$ than the basis has entries, @f$ \boldsymbol{M} @f$ is singular;
   * the completeness order is then reduced for this evaluation point (getCorrectedCompletenessOrder()).
   *
   * **Limit on the condition of the moment matrix (optional).** Enough kernels may cover @f$ \boldsymbol{x} @f$ and
   * still leave @f$ \boldsymbol{M} @f$ nearly singular: at a point on one side of all covering kernels, each at the
   * edge of its support (e.g. the corner of a strongly stretched particle reaching into a void), the corrected
   * kernels extrapolate and @f$ \Psi_A @f$ and its gradients become arbitrarily large. With a limit
   * @f$ \kappa_\mathrm{max} @f$ on the condition of @f$ \boldsymbol{M} @f$, the completeness order of such a point is
   * lowered one order at a time until the condition is below the limit, down to the order @f$ 0 @f$ (Shepard
   * functions, @f$ \Psi_A = \phi_A / \sum_B \phi_B @f$, bounded by 1), which is always accepted. The condition is that
   * of the diagonally scaled matrix @f$ \boldsymbol{D} \boldsymbol{M} \boldsymbol{D} @f$,
   * @f$ D_{kk} = M_{kk}^{-1/2} @f$ (conditionOfMomentMatrix()): the scaling makes it independent of the unit of length
   * and of a constant factor of the kernels. Without a limit (the default, @f$ \kappa_\mathrm{max} = \infty @f$) the
   * condition is not computed; with a limit, a point below it is evaluated exactly as without one, and a moment
   * matrix of order @f$ n \ge 1 @f$ that is singular is evaluated with a lower order as well (only a point that no
   * kernel covers still throws). The lower order costs the point its polynomial completeness: the shape functions
   * jump where the condition crosses the limit, in space and in time, and a particle integrated with such a point no
   * longer passes the patch test of that order.
   *
   * See ReproducingKernelApproximationImplicit for the variant with implicit gradients.
   */
  class ReproducingKernelApproximation : public Approximation {

  public:
    /// The largest spatial dimension, that of the kernels.
    static constexpr int maxDim = KernelFunction::maxDim;
    /// The largest monomial basis: the complete cubic basis in 3D, @f$ \binom{3 + 3}{3} @f$.
    static constexpr int maxSizeOfBasis = 20;

    /// The monomial basis @f$ \boldsymbol{H} @f$ (and vectors of its size), of a size up to maxSizeOfBasis.
    using BasisVector = Eigen::Matrix< double, Eigen::Dynamic, 1, Eigen::ColMajor, maxSizeOfBasis, 1 >;
    /// The gradient @f$ \partial \boldsymbol{H} / \partial \boldsymbol{x} @f$, one column per dimension.
    using BasisGradient = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, maxSizeOfBasis, maxDim >;
    /// The moment matrix @f$ \boldsymbol{M} @f$ (and matrices of its size).
    using MomentMatrix = Eigen::
      Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, maxSizeOfBasis, maxSizeOfBasis >;
    /// The gradients @f$ \boldsymbol{M}_{,i} @f$ of the moment matrix, entry @f$ i @f$.
    using MomentMatrixGradient = std::array< MomentMatrix, maxDim >;
    /// A point, or a difference of points, of up to maxDim coordinates.
    using Coordinates = Eigen::Matrix< double, Eigen::Dynamic, 1, Eigen::ColMajor, maxDim, 1 >;

  protected:
    int _dim;                      ///< spatial dimension @f$ d @f$
    int _desiredCompletenessOrder; ///< requested completeness order @f$ n @f$ (may be reduced per point)
    /// The limit @f$ \kappa_\mathrm{max} @f$ on the condition of the scaled moment matrix, above which a point is
    /// evaluated with the completeness order 0; infinite for no limit.
    double _momentMatrixConditionLimit;
    /// The monomial basis of every order up to the desired one, entry @f$ n @f$ (an order may be reduced per point).
    std::vector< Math::MonomialBasis > _basisOfOrder;

    /**
     * @brief Evaluate every candidate kernel at a point, once: its value @f$ \phi_A(\boldsymbol{x}) @f$ and, if
     * asked for, its gradient. A kernel with @f$ |\phi_A(\boldsymbol{x})| \le 10^{-14} @f$ does not cover the
     * point: its value (and gradient) is set to zero, and it takes no further part.
     * @param[in]  coord           Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] phi             The kernel values, length @f$ n_\mathrm{c} @f$.
     * @param[out] dPhi_dx         The kernel gradients, column-major @f$ d \times n_\mathrm{c} @f$; nullptr if not
     *                             needed.
     * @return The number of kernels that cover the point.
     */
    int evaluateKernelFunctions( const double*                               coord,
                                 const std::vector< const KernelFunction* >& kernelFunctions,
                                 double*                                     phi,
                                 double*                                     dPhi_dx ) const;

    /**
     * @brief The monomial basis @f$ \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_A) @f$ of every covering kernel at
     * its shifted coordinates, and, if asked for, its gradient (Math::MonomialBasis): evaluated once per kernel and
     * point, for the moment matrix and the shape functions alike. Since @f$ \boldsymbol{x}_A @f$ is fixed, the gradient
     * with respect to the shifted coordinates is the one with respect to @f$ \boldsymbol{x} @f$.
     * @param[in]  coord             Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in]  kernelFunctions   The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[in]  phi               Their values, from evaluateKernelFunctions(); a kernel with zero value does not
     *                               cover the point and is skipped (its columns are left as they are).
     * @param[in]  completenessOrder Completeness order @f$ n @f$.
     * @param[out] H                 The bases, column @f$ A @f$ is @f$ \boldsymbol{H}(\boldsymbol{x} -
     *                               \boldsymbol{x}_A) @f$; resized to the size of the basis times @f$ n_\mathrm{c} @f$.
     * @param[out] dH_dx             Their gradients, column @f$ d A + i @f$ is @f$ \boldsymbol{H}_{,i} @f$ of kernel
     *                               @f$ A @f$; resized to the size of the basis times @f$ d\, n_\mathrm{c} @f$; nullptr
     *                               if not needed.
     */
    void evaluateBasisOfKernels( const double*                               coord,
                                 const std::vector< const KernelFunction* >& kernelFunctions,
                                 const double*                               phi,
                                 int                                         completenessOrder,
                                 Eigen::MatrixXd&                            H,
                                 Eigen::MatrixXd*                            dH_dx ) const;

    /**
     * @brief Assemble the moment matrix @f$ \boldsymbol{M}(\boldsymbol{x}) = \sum_B \boldsymbol{H}(\boldsymbol{x} -
     * \boldsymbol{x}_B) \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_B) \phi_B(\boldsymbol{x}) @f$ over the
     * covering kernels and, if asked for, its gradients @f$ \boldsymbol{M}_{,i} = \sum_B ( \boldsymbol{H}_{,i}
     * \boldsymbol{H}^T + \boldsymbol{H} \boldsymbol{H}_{,i}^T ) \phi_B + \boldsymbol{H} \boldsymbol{H}^T
     * \phi_{B,i} @f$: the one assembly of the moment matrix.
     * @param[in]  H       The basis of every candidate kernel, column @f$ B @f$ (evaluateBasisOfKernels()).
     * @param[in]  dH_dx   Its gradients, columns @f$ d B + i @f$; nullptr if @p dM_dx is.
     * @param[in]  phi     The kernel values, from evaluateKernelFunctions(); zero for a kernel that does not cover
     *                     the point, which is skipped.
     * @param[in]  dPhi_dx The kernel gradients, from evaluateKernelFunctions(); nullptr if @p dM_dx is.
     * @param[out] M       The symmetric moment matrix, of the size of the basis squared.
     * @param[out] dM_dx   Its @f$ d @f$ gradients; nullptr if not needed.
     */
    void assembleMomentMatrix( const Eigen::MatrixXd& H,
                               const Eigen::MatrixXd* dH_dx,
                               const double*          phi,
                               const double*          dPhi_dx,
                               MomentMatrix&          M,
                               MomentMatrixGradient*  dM_dx ) const;

    /**
     * @brief Factorize the moment matrix @f$ \boldsymbol{M}(\boldsymbol{x}) @f$ and check that it is regular.
     * @details checkNonSingularity() only compares the number of covering kernels with the size of the basis, which
     *          is necessary but not sufficient: for nodes in a degenerate arrangement (e.g. collinear nodes in 2D at
     *          order 1), or for a point covered by no kernel at all, @f$ \boldsymbol{M} @f$ is singular, and solving
     *          it would silently give shape functions without partition of unity.
     * @param[in] M                 Moment matrix.
     * @param[in] coord             Evaluation point (for the error message).
     * @param[in] nCoveringKernels  Number of kernels covering the point (for the error message).
     * @return The column-pivoting Householder QR factorization of @p M.
     * @throws std::runtime_error if @p M is rank deficient. A common cause in 3D is a single layer of particles
     *         through the thickness: all kernel centres are coplanar, so a linear basis cannot be reproduced; use at
     *         least two layers of kernels, larger supports, or a lower completeness order.
     */
    Eigen::ColPivHouseholderQR< MomentMatrix > factorizeMomentMatrix( const MomentMatrix& M,
                                                                      const double*       coord,
                                                                      int                 nCoveringKernels ) const;

    /**
     * @brief Choose the completeness order of a point, and assemble its moment matrix of that order: the desired
     * order, reduced for few covering kernels (getCorrectedCompletenessOrder()), and lowered one order at a time while
     * the condition of the moment matrix exceeds the limit @f$ \kappa_\mathrm{max} @f$ (conditionOfMomentMatrix()),
     * down to the order 0, which is always accepted; the one choice of the order for both approximations. Without a
     * limit, or below it, this is evaluateBasisOfKernels() and assembleMomentMatrix() at the order for the number of
     * covering kernels.
     * @param[in]  coord           Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in]  kernelFunctions The candidate kernel functions.
     * @param[in]  phi             Their values, from evaluateKernelFunctions().
     * @param[in]  dPhi_dx         Their gradients, from evaluateKernelFunctions(); nullptr if @p dM_dx is.
     * @param[in]  nCovering       The number of kernels that cover the point, from evaluateKernelFunctions().
     * @param[out] H               The bases of the kernels at the chosen order (evaluateBasisOfKernels()).
     * @param[out] dH_dx           Their gradients; nullptr if not needed.
     * @param[out] M               The moment matrix of the chosen order.
     * @param[out] dM_dx           Its @f$ d @f$ gradients; nullptr if not needed.
     * @return The completeness order the point is evaluated with.
     */
    int chooseCompletenessOrderAndAssembleMomentMatrix( const double*                               coord,
                                                        const std::vector< const KernelFunction* >& kernelFunctions,
                                                        const double*                               phi,
                                                        const double*                               dPhi_dx,
                                                        int                                         nCovering,
                                                        Eigen::MatrixXd&                            H,
                                                        Eigen::MatrixXd*                            dH_dx,
                                                        MomentMatrix&                               M,
                                                        MomentMatrixGradient*                       dM_dx ) const;

    /**
     * @brief The condition of the moment matrix, scaled to a unit diagonal: the estimate of the condition number in
     * the 1-norm, @f$ \kappa_1(\boldsymbol{D} \boldsymbol{M} \boldsymbol{D}) @f$ with @f$ D_{kk} = M_{kk}^{-1/2}
     * @f$, from a Cholesky (LDLT) factorization of the scaled matrix and its reciprocal condition estimate.
     *
     * The monomials of @f$ \boldsymbol{H} @f$ carry the unit of length to their degree, so the condition of
     * @f$ \boldsymbol{M} @f$ itself would depend on the unit; scaling to a unit diagonal removes this, and a
     * constant factor of the kernels. In the 2-norm, the scaling to a unit diagonal gives a condition within a factor
     * of the basis size @f$ \binom{d+n}{d} @f$ of the best diagonal scaling (van der Sluis); the 1-norm condition
     * differs from the 2-norm one by at most that factor as well. For the small, symmetric moment matrices the
     * factorization costs less than the QR factorization of the evaluation.
     * @param[in] M The moment matrix.
     * @return @f$ \kappa_1 \ge 1 @f$; infinite for a singular matrix (a pivot of the scaled matrix at the rounding
     *         level) or a zero entry on its diagonal.
     */
    static double conditionOfMomentMatrix( const MomentMatrix& M );

    /**
     * @brief Check if the completeness order leads to a non-singular equation system.
     *
     * Checks the necessary condition that at least as many kernels cover the point as the basis has entries,
     * @f$ n_\mathrm{nodes} \ge \binom{d + n}{d} @f$, the size of the basis. It is not sufficient, e.g., for nodes in a
     * degenerate (collinear, coplanar) arrangement. Always true for @f$ n \le 0 @f$.
     *
     * @param[in] nNodes            Number of kernels covering the evaluation point.
     * @param[in] completenessOrder Completeness order @f$ n @f$.
     * @return True if the condition is met.
     */
    bool checkNonSingularity( int nNodes, int completenessOrder ) const
    {
      return completenessOrder <= 0 || nNodes >= _basisOfOrder[completenessOrder].size();
    }

    /**
     * @brief Compute the corrected completeness order.
     *
     * Depending on the actual number of covering kernel functions, the desired completeness order might not be
     * possible. Starting from the desired order, the order is decreased until checkNonSingularity() holds.
     *
     * @param[in] nNodes Number of kernels covering the evaluation point.
     * @return The highest order @f$ \le @f$ the desired one for which checkNonSingularity() holds.
     */
    int getCorrectedCompletenessOrder( int nNodes ) const
    {
      int correctedCompletenessOrder = _desiredCompletenessOrder;
      while ( !checkNonSingularity( nNodes, correctedCompletenessOrder ) )
        correctedCompletenessOrder--;
      return correctedCompletenessOrder;
    }

  public:
    /**
     * @brief Construct the approximation.
     * @param[in] dim                        Spatial dimension @f$ d \le @f$ maxDim.
     * @param[in] completenessOrder          Desired completeness order @f$ n \ge 0 @f$ (order of the reproduced
     *                                       polynomials).
     * @param[in] momentMatrixConditionLimit The limit @f$ \kappa_\mathrm{max} \ge 1 @f$ on the condition of the
     *                                       scaled moment matrix, above which the completeness order of a point
     *                                       is lowered; infinite (the default) for no limit.
     * @throws std::invalid_argument if the order is negative, the dimension or the basis exceeds maxDim or
     * maxSizeOfBasis, or the limit is below 1 or not a number.
     */
    ReproducingKernelApproximation( int    dim,
                                    int    completenessOrder,
                                    double momentMatrixConditionLimit = std::numeric_limits< double >::infinity() );

    /// @brief Default destructor.
    virtual ~ReproducingKernelApproximation() = default;

    /**
     * @brief Compute the RK shape function values
     * @f$ \Psi_A = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A) \boldsymbol{b} \phi_A @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     * @throws std::runtime_error if the moment matrix is singular (see factorizeMomentMatrix()).
     */
    virtual void computeShapeFunctions( const double*                               coord,
                                        const std::vector< const KernelFunction* >& kernelFunctions,
                                        double* shapeFunctionValues ) const override;

    /**
     * @brief Compute the RK shape function values and their explicit (exact) gradients @f$ \Psi_{A,i} @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     * @param[out] shapeFunctionValueGradients_ColMajor Gradients, column-major @f$ d \times n_\mathrm{c} @f$
     *                                 (entry @f$ i + d A @f$ is @f$ \Psi_{A,i} @f$); zero for non-covering candidates.
     * @throws std::runtime_error if the moment matrix is singular (see factorizeMomentMatrix()).
     */
    virtual void computeShapeFunctionsAndGradients( const double*                               coord,
                                                    const std::vector< const KernelFunction* >& kernelFunctions,
                                                    double*                                     shapeFunctionValues,
                                                    double* shapeFunctionValueGradients_ColMajor ) const override;

  private:
    /**
     * @brief The evaluation behind computeShapeFunctions() and computeShapeFunctionsAndGradients(): the values,
     * and the gradients if @p dPsi_dx is given.
     * @param[in]  coord           Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in]  kernelFunctions The candidate kernel functions.
     * @param[out] Psi             Shape function values.
     * @param[out] dPsi_dx         Their gradients, column-major @f$ d \times n_\mathrm{c} @f$; nullptr for the values
     *                             alone.
     */
    void evaluateShapeFunctions( const double*                               coord,
                                 const std::vector< const KernelFunction* >& kernelFunctions,
                                 double*                                     Psi,
                                 double*                                     dPsi_dx ) const;
  };

} // namespace Marmot::Discretizations::Meshfree
