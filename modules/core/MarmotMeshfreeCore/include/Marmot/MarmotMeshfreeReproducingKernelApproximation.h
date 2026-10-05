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

#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include <Eigen/Core>
#include <Eigen/QR>
#include <cmath>
#include <vector>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::MarmotMeshfreeReproducingKernelApproximation
   * @brief Reproducing kernel (RK) approximation with explicit (direct) shape function gradients.
   *
   * The shape function of node @f$ A @f$ is the kernel @f$ \phi_A @f$ (MarmotMeshfreeKernelFunction) multiplied
   * by a correction function,
   * @f[
   *   \Psi_A(\boldsymbol{x}) = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A)\, \boldsymbol{b}(\boldsymbol{x})\,
   *   \phi_A(\boldsymbol{x}),
   * @f]
   * where @f$ \boldsymbol{H} @f$ is the vector of all monomials up to the completeness order @f$ n @f$
   * (Math::computeMonomialBasis(), first entry @f$ 1 @f$) evaluated at the shifted coordinates
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
   * If fewer kernels cover @f$ \boldsymbol{x} @f$ than the basis has entries, @f$ \boldsymbol{M} @f$ is singular;
   * the completeness order is then reduced for this evaluation point (getCorrectedCompletenessOrder()).
   *
   * See MarmotMeshfreeReproducingKernelApproximationImplicit for the variant with implicit gradients.
   */
  class MarmotMeshfreeReproducingKernelApproximation : public MarmotMeshfreeApproximation {

  protected:
    int _dim;                      ///< spatial dimension @f$ d @f$
    int _desiredCompletenessOrder; ///< requested completeness order @f$ n @f$ (may be reduced per point)

    // int static computeHRecursively( int                    completenessOrder,
    //                                 const Eigen::VectorXd& x_minus_center,
    //                                 Eigen::VectorXd&       res,
    //                                 int                    idx,
    //                                 int                    dim );

    // int static computeHGradientRecursively( int                    completenessOrder,
    //                                         const Eigen::VectorXd& x_minus_center,
    //                                         Eigen::MatrixXd&       res,
    //                                         int                    idx,
    //                                         int                    dim );

    // static int computeSizeHVector( int completenessOrder, int dim );

    /**
     * @brief Compute the monomial basis vector @f$ \boldsymbol{H} @f$ at shifted coordinates.
     *
     * Wraps Math::computeMonomialBasis(); callers pass @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$.
     *
     * @param[in] globalCoord            Shifted coordinates @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$ at which the
     *                                   basis is evaluated (despite the name, not the global coordinates).
     * @param[in] coveringShapeFunctions Covering kernel functions (unused).
     * @param[in] completenessOrder      Completeness order @f$ n @f$.
     * @return @f$ \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_A) @f$, of size
     *         @f$ \binom{d + n}{d} @f$.
     */
    static Eigen::VectorXd computeHVector(
      const Eigen::VectorXd&                                    globalCoord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
      int                                                       completenessOrder );

    /**
     * @brief Compute the gradient @f$ \partial \boldsymbol{H} / \partial x_i @f$ of the monomial basis vector.
     *
     * Wraps Math::computeMonomialBasisGradient(); callers pass @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$. Since
     * @f$ \boldsymbol{x}_A @f$ is fixed, the derivative with respect to the shifted coordinates equals that with
     * respect to @f$ \boldsymbol{x} @f$.
     *
     * @param[in] globalCoord            Shifted coordinates @f$ \boldsymbol{x} - \boldsymbol{x}_A @f$ (despite the
     *                                   name, not the global coordinates).
     * @param[in] coveringShapeFunctions Covering kernel functions (unused).
     * @param[in] completenessOrder      Completeness order @f$ n @f$.
     * @return Matrix of size @f$ \binom{d + n}{d} \times d @f$; column @f$ i @f$ holds
     *         @f$ \partial \boldsymbol{H} / \partial x_i @f$.
     */
    static Eigen::MatrixXd computeHVectorGradient(
      const Eigen::VectorXd&                                    globalCoord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
      int                                                       completenessOrder );

    /**
     * @brief Compute the moment matrix
     * @f$ \boldsymbol{M}(\boldsymbol{x}) = \sum_B \boldsymbol{H}(\boldsymbol{x} - \boldsymbol{x}_B)
     * \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_B) \phi_B(\boldsymbol{x}) @f$.
     * @param[in] globalCoord            Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in] coveringShapeFunctions Kernel functions @f$ \phi_B @f$ summed over (kernels that do not cover
     *                                   @f$ \boldsymbol{x} @f$ contribute zero).
     * @param[in] completenessOrder      Completeness order @f$ n @f$.
     * @return The symmetric moment matrix, of size @f$ \binom{d + n}{d} \times \binom{d + n}{d} @f$.
     */
    static Eigen::MatrixXd computeMMatrix(
      const Eigen::VectorXd&                                    globalCoord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
      int                                                       completenessOrder );

    /**
     * @brief Compute the moment matrix @f$ \boldsymbol{M} @f$ and its gradients @f$ \boldsymbol{M}_{,i} @f$.
     *
     * @f$ \boldsymbol{M}_{,i} = \sum_B ( \boldsymbol{H}_{,i} \boldsymbol{H}^T + \boldsymbol{H}
     * \boldsymbol{H}_{,i}^T ) \phi_B + \boldsymbol{H} \boldsymbol{H}^T \phi_{B,i} @f$, with
     * @f$ \boldsymbol{H} @f$ evaluated at @f$ \boldsymbol{x} - \boldsymbol{x}_B @f$.
     *
     * @param[in] globalCoord            Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in] coveringShapeFunctions Kernel functions @f$ \phi_B @f$ summed over.
     * @param[in] completenessOrder      Completeness order @f$ n @f$.
     * @return Pair of @f$ \boldsymbol{M} @f$ and the @f$ d @f$ matrices @f$ \boldsymbol{M}_{,i} @f$.
     */
    static std::pair< Eigen::MatrixXd, std::vector< Eigen::MatrixXd > > computeMMatrixAndGradient(
      const Eigen::VectorXd&                                    globalCoord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& coveringShapeFunctions,
      int                                                       completenessOrder );

    /**
     * @brief Get @f$ \boldsymbol{H}_0 = \boldsymbol{H}(\boldsymbol{0}) = [1, 0, \dots, 0]^T @f$.
     * @param[in] sizeHVector Size of the monomial basis vector.
     * @return The first unit vector of the given size.
     */
    static Eigen::VectorXd H0Vector( int sizeHVector );

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
     * @throws std::runtime_error if @p M is rank deficient.
     */
    Eigen::ColPivHouseholderQR< Eigen::MatrixXd > factorizeMomentMatrix( const Eigen::MatrixXd& M,
                                                                         const double*          coord,
                                                                         int nCoveringKernels ) const;

    /**
     * @brief Compute the factorial @f$ n! @f$ recursively.
     * @param[in] n Non-negative integer.
     * @return @f$ n! @f$.
     */
    int factorial( int n ) const
    {
      if ( n == 0 )
        return 1;
      else
        return n * factorial( n - 1 );
    }

    /**
     * @brief Check if the completeness order leads to a non-singular equation system.
     *
     * Checks the necessary condition that at least as many kernels cover the point as the basis has entries,
     * @f$ n_\mathrm{nodes} \ge \binom{d + n}{d} @f$. It is not sufficient, e.g., for nodes in a degenerate
     * (collinear, coplanar) arrangement. Always true for @f$ n \le 0 @f$.
     *
     * @param[in] nNodes            Number of kernels covering the evaluation point.
     * @param[in] completenessOrder Completeness order @f$ n @f$.
     * @return True if the condition is met.
     */
    bool checkNonSingularity( int nNodes, int completenessOrder ) const
    {
      if ( completenessOrder > 0 ) {
        return nNodes >= factorial( _dim + completenessOrder ) / ( factorial( _dim ) * factorial( completenessOrder ) );
      }
      else
        return true;
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

    /**
     * @brief Find the kernels covering a point.
     * @param[in] coord           Evaluation point @f$ \boldsymbol{x} @f$.
     * @param[in] kernelFunctions Candidate kernel functions.
     * @return Indices (into @p kernelFunctions) of the kernels with @f$ |\phi_A(\boldsymbol{x})| > 10^{-14} @f$.
     */
    const std::vector< int > findCoveringKernelFunctionIndices(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions ) const;

  public:
    /**
     * @brief Construct the approximation.
     * @param[in] dim               Spatial dimension @f$ d @f$.
     * @param[in] completenessOrder Desired completeness order @f$ n @f$ (order of the reproduced polynomials).
     */
    MarmotMeshfreeReproducingKernelApproximation( int dim, int completenessOrder );

    /// @brief Default destructor.
    virtual ~MarmotMeshfreeReproducingKernelApproximation() = default;

    /**
     * @brief Compute the RK shape function values
     * @f$ \Psi_A = \boldsymbol{H}^T(\boldsymbol{x} - \boldsymbol{x}_A) \boldsymbol{b} \phi_A @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     */
    virtual void computeShapeFunctions( const double*                                             coord,
                                        const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
                                        double* shapeFunctionValues ) const override;

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
     * @brief Compute the RK shape function values and their explicit (exact) gradients @f$ \Psi_{A,i} @f$.
     * @param[in]  coord               Evaluation point @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     * @param[in]  kernelFunctions     The @f$ n_\mathrm{c} @f$ candidate kernel functions.
     * @param[out] shapeFunctionValues Shape function values (length @f$ n_\mathrm{c} @f$); zero for candidates that
     *                                 do not cover @f$ \boldsymbol{x} @f$.
     * @param[out] shapeFunctionValueGradients_ColMajor Gradients, column-major @f$ d \times n_\mathrm{c} @f$
     *                                 (entry @f$ i + d A @f$ is @f$ \Psi_{A,i} @f$); zero for non-covering candidates.
     * @throws std::runtime_error if the basis is empty (negative completeness order).
     */
    virtual void computeShapeFunctionsAndGradients(
      const double*                                             coord,
      const std::vector< const MarmotMeshfreeKernelFunction* >& kernelFunctions,
      double*                                                   shapeFunctionValues,
      double*                                                   shapeFunctionValueGradients_ColMajor ) const override;
  };

} // namespace Marmot::Meshfree
