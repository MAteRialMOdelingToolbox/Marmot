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
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include <cmath>
#include <functional>
#include <tuple>

namespace Marmot::Math {

  /** @brief Value and first/second derivative w.r.t. a symmetric 3x3 tensor \f$\boldsymbol A\f$ of an
   * isotropic "sum-of-eigenvalue-function" \f$g(\boldsymbol A) = \sum_{i=1}^3 h(\lambda_i)\f$, correct and
   * numerically robust EVEN AT repeated eigenvalues of \f$\boldsymbol A\f$ (including the fully isotropic
   * case \f$\boldsymbol A=c\boldsymbol I\f$).
   *
   * @details This is the closed-form building block needed by any hyperelastic potential built from a sum
   * of a function of \f$\boldsymbol A\f$'s principal values (e.g. an Ogden-type potential
   * \f$\Psi=\sum_p \frac{\mu_p}{\alpha_p}\sum_i\left(\bar\lambda_i^{\alpha_p}-1\right)\f$, which is exactly
   * of this form after the substitution \f$\boldsymbol{\bar A} = \det(\boldsymbol A)^{-1/3}\boldsymbol A\f$
   * and \f$h(\mu) = \sum_p \frac{\mu_p}{\alpha_p}\mu^{\alpha_p/2}\f$, since \f$\bar\lambda_i^2 =
   * \mathrm{eig}_i(\boldsymbol{\bar A})\f$).
   *
   * @par Why this exists (background)
   * A naive way to differentiate such a potential is to run a generic eigen-decomposition (e.g.
   * #computeEigenSystemJacobi) with an automatic-differentiation scalar type and let the AD machinery
   * propagate through the decomposition's own arithmetic. This FAILS SILENTLY at repeated eigenvalues
   * (most notably \f$\boldsymbol A=\boldsymbol I\f$, i.e. every reference/undeformed configuration): Jacobi's
   * convergence and per-pair rotation gates are keyed on the PRIMAL magnitude of the off-diagonal entries
   * (see #computeEigenSystemJacobi's use of #makeReal), so any derivative content carried by a dual or
   * complex-step number whose PRIMAL value is (numerically) zero -- exactly the case for every off-diagonal
   * entry at a repeated eigenvalue -- is silently discarded rather than propagated into the eigenvalues,
   * even though the true scalar \f$g(\boldsymbol A)\f$, being a permutation-symmetric function of the full
   * eigenvalue set, IS mathematically smooth there. Relaxing the convergence gates naively does not fix
   * this either: it was checked and found to introduce literal primal-level \f$0/0\f$ divisions (NaN) for
   * genuinely uncoupled off-diagonal pairs.
   *
   * @par The fix
   * This function instead computes the derivatives via the classical closed-form "isotropic tensor
   * function" / Daleckii-Krein divided-difference formulas (see e.g. M. Itskov, "Tensor Algebra and Tensor
   * Analysis for Engineers", or C. Miehe (1998), "Comparison of two algorithms for the computation of
   * fourth-order isotropic tensor functions", Comput. Struct. 66(1):37-43, or de Souza Neto et al.,
   * "Computational Methods for Plasticity", Box 12.1), operating ENTIRELY on the PRIMAL (double) eigenvalue
   * decomposition -- no automatic differentiation through the decomposition is ever needed. With
   * eigenvalues \f$\lambda_i\f$ and eigenprojectors \f$\boldsymbol E_i=\boldsymbol n_i\otimes\boldsymbol
   * n_i\f$ (\f$\boldsymbol A=\sum_i\lambda_i\boldsymbol E_i\f$):
   * \f[
   *   \frac{\partial g}{\partial \boldsymbol A} = \sum_{i=1}^3 h'(\lambda_i)\,\boldsymbol E_i,
   *   \qquad
   *   \frac{\partial^2 g}{\partial \boldsymbol A\,\partial \boldsymbol A} = \sum_{i,j=1}^3
   *   h'[\lambda_i,\lambda_j]\,\left(\boldsymbol E_i \odot \boldsymbol E_j\right)
   * \f]
   * where \f$(\boldsymbol X\odot\boldsymbol Y)_{abcd}=\frac12(X_{ac}Y_{bd}+X_{ad}Y_{bc})\f$ is the minor-
   * symmetric "square" tensor product, and \f$h'[\lambda_i,\lambda_j]\f$ is the first divided difference of
   * \f$h'\f$,
   * \f[
   *   h'[\lambda_i,\lambda_j] = \begin{cases} \dfrac{h'(\lambda_i)-h'(\lambda_j)}{\lambda_i-\lambda_j} &
   *   \lambda_i\neq\lambda_j \\[4pt] h''(\lambda_i) & \lambda_i=\lambda_j\end{cases}
   * \f]
   * -- the \f$\lambda_i=\lambda_j\f$ branch is exactly the removable-singularity (L'Hopital) limit of the
   * generic formula, so \f$h'[\cdot,\cdot]\f$ (and therefore the whole second derivative) is smooth across
   * the repeated-eigenvalue transition BY CONSTRUCTION -- unlike differentiating the decomposition itself,
   * no information is ever discarded. Crucially, the first derivative formula above uses no divided
   * differences at all and is valid at ANY eigenvalue multiplicity without special-casing, since summing
   * \f$h'(\lambda_i)\boldsymbol E_i\f$ over a repeated cluster is basis-independent (equal to
   * \f$h'(\lambda)\f$ times the projector onto that eigenspace, regardless of which orthonormal eigenbasis
   * was chosen within it).
   *
   * @par Validation
   * Verified standalone (independent of any specific hyperelastic potential) against: (1) the exactly known
   * closed forms for \f$h(\lambda)=\lambda\f$ (\f$g=\mathrm{tr}\boldsymbol A\f$) and \f$h(\lambda)=\lambda^2\f$
   * (\f$g=\mathrm{tr}(\boldsymbol A^2)\f$); and (2) raw finite differences of \f$g\f$ and of this function's
   * own first derivative, for both polynomial and non-integer-exponent \f$h\f$ (mimicking a real Ogden
   * exponent), at distinct, partially repeated, and fully repeated (\f$\boldsymbol A=c\boldsymbol I\f$)
   * eigenvalue configurations.
   *
   * @param A Symmetric tensor at which to evaluate \f$g\f$ and its derivatives.
   * @param h Scalar function applied to each eigenvalue.
   * @param hPrime First derivative of @p h.
   * @param hDoublePrime Second derivative of @p h.
   * @param degenerateEigenvalueTol Absolute eigenvalue-difference threshold below which the L'Hopital limit
   * \f$h''(\lambda_i)\f$ is used in place of the (numerically ill-conditioned) divided-difference quotient.
   * @return Tuple of \f$\{g,\ \partial g/\partial\boldsymbol A,\ \partial^2 g/\partial\boldsymbol
   * A\partial\boldsymbol A\}\f$.
   */
  inline std::tuple< double, FastorStandardTensors::Tensor33d, FastorStandardTensors::Tensor3333d > sumOfEigenvaluesFunctionAndDerivatives(
    const FastorStandardTensors::Tensor33d&  A,
    const std::function< double( double ) >& h,
    const std::function< double( double ) >& hPrime,
    const std::function< double( double ) >& hDoublePrime,
    double                                   degenerateEigenvalueTol = 1e-6 )
  {
    using namespace FastorStandardTensors;
    using namespace Fastor;
    using namespace FastorIndices;

    const auto [lambda, N] = computeEigenSystemJacobi< double >( A );

    Tensor33d E[3];
    for ( int i = 0; i < 3; ++i ) {
      Tensor3d ni;
      for ( int k = 0; k < 3; ++k )
        ni( k ) = N( k, i );
      E[i] = outer( ni, ni );
    }

    auto firstDividedDifferenceOfHPrime = [&]( double a, double b ) {
      if ( std::abs( a - b ) < degenerateEigenvalueTol )
        return hDoublePrime( 0.5 * ( a + b ) );
      return ( hPrime( a ) - hPrime( b ) ) / ( a - b );
    };

    auto squareProduct = []( const Tensor33d& X, const Tensor33d& Y ) {
      return Tensor3333t< double >(
        evaluate( 0.5 * ( einsum< ik, jl, to_ijkl >( X, Y ) + einsum< il, jk, to_ijkl >( X, Y ) ) ) );
    };

    double g = 0.0;
    for ( int i = 0; i < 3; ++i )
      g += h( lambda( i ) );

    Tensor33d dg_dA( 0.0 );
    for ( int i = 0; i < 3; ++i )
      dg_dA += hPrime( lambda( i ) ) * E[i];

    Tensor3333d d2g_dAdA( 0.0 );
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        d2g_dAdA += firstDividedDifferenceOfHPrime( lambda( i ), lambda( j ) ) * squareProduct( E[i], E[j] );

    return { g, dg_dA, d2g_dAdA };
  }

  /** @brief As #sumOfEigenvaluesFunctionAndDerivatives, but additionally returns the third derivative
   * w.r.t. \f$\boldsymbol A\f$.
   *
   * @details The third derivative is obtained by central-differencing the (already closed-form and
   * repeated-eigenvalue-safe) SECOND derivative #sumOfEigenvaluesFunctionAndDerivatives itself, rather than
   * by hand-deriving the considerably more involved general triple-eigenprojector closed form: since the
   * second derivative has no poles left (the divided-difference construction already removed them),
   * differentiating it with an ordinary real perturbation is a ordinary, well-behaved numerical operation
   * -- a generic small real perturbation of \f$\boldsymbol A\f$ does not reintroduce the exact-zero-primal
   * situation that broke the naive AD-through-eigendecomposition approach in the first place.
   *
   * Validated standalone against an INDEPENDENT raw directional triple finite difference of \f$g\f$ itself
   * (the standard 8-corner stencil, bypassing this function's own first/second derivative machinery
   * entirely), for direction triples deliberately including pure off-diagonal perturbations at a fully
   * repeated-eigenvalue \f$\boldsymbol A\f$ -- i.e. exactly the configuration that silently broke under the
   * old approach. Relative agreement was consistently within the expected finite-difference truncation
   * error of the step size used here (better than \f$10^{-5}\f$ for @p thirdDerivativeFDStep \f$=10^{-5}\f$).
   *
   * @param A Symmetric tensor at which to evaluate \f$g\f$ and its derivatives.
   * @param h Scalar function applied to each eigenvalue.
   * @param hPrime First derivative of @p h.
   * @param hDoublePrime Second derivative of @p h.
   * @param degenerateEigenvalueTol Forwarded to #sumOfEigenvaluesFunctionAndDerivatives.
   * @param thirdDerivativeFDStep Central-difference step size used for the third derivative. Smaller is not
   * necessarily better (finite-difference truncation vs. floating-point cancellation trade-off); the
   * default was chosen empirically and matched the validation above.
   * @return Tuple of \f$\{g,\ \partial g/\partial\boldsymbol A,\ \partial^2 g/\partial\boldsymbol
   * A\partial\boldsymbol A,\ \partial^3 g/\partial\boldsymbol A\partial\boldsymbol A\partial\boldsymbol
   * A\}\f$.
   */
  inline std::tuple< double,
                     FastorStandardTensors::Tensor33d,
                     FastorStandardTensors::Tensor3333d,
                     FastorStandardTensors::Tensor333333d >
  sumOfEigenvaluesFunctionAndDerivativesUpToThird( const FastorStandardTensors::Tensor33d&  A,
                                                   const std::function< double( double ) >& h,
                                                   const std::function< double( double ) >& hPrime,
                                                   const std::function< double( double ) >& hDoublePrime,
                                                   double degenerateEigenvalueTol = 1e-6,
                                                   double thirdDerivativeFDStep   = 1e-5 )
  {
    using namespace FastorStandardTensors;

    const auto [g, dg_dA, d2g_dAdA] = sumOfEigenvaluesFunctionAndDerivatives( A,
                                                                              h,
                                                                              hPrime,
                                                                              hDoublePrime,
                                                                              degenerateEigenvalueTol );

    auto symmetricallyPerturbed = [&]( int m, int n, double delta ) {
      Tensor33d Ap = A;
      Ap( m, n ) += delta;
      if ( m != n )
        Ap( n, m ) += delta;
      return Ap;
    };

    Tensor333333d d3g_dAdAdA( 0.0 );
    for ( int m = 0; m < 3; ++m ) {
      for ( int n = 0; n < 3; ++n ) {
        const double delta = ( m == n ) ? thirdDerivativeFDStep : thirdDerivativeFDStep / 2.0;

        const auto [gp, dgp, d2p] = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( m, n, delta ),
                                                                            h,
                                                                            hPrime,
                                                                            hDoublePrime,
                                                                            degenerateEigenvalueTol );
        const auto [gm, dgm, d2m] = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( m, n, -delta ),
                                                                            h,
                                                                            hPrime,
                                                                            hDoublePrime,
                                                                            degenerateEigenvalueTol );

        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            for ( int k = 0; k < 3; ++k )
              for ( int l = 0; l < 3; ++l )
                d3g_dAdAdA( i, j, k, l, m, n ) = ( d2p( i, j, k, l ) - d2m( i, j, k, l ) ) /
                                                 ( 2.0 * thirdDerivativeFDStep );
      }
    }

    return { g, dg_dA, d2g_dAdA, d3g_dAdAdA };
  }

} // namespace Marmot::Math
