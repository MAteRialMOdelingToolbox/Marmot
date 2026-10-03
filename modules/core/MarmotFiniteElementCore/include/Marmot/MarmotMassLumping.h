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
#include <algorithm>
#include <cmath>
#include <limits>

namespace Marmot::FiniteElement::MassLumping {

  /**
   * @brief The blend weight of the manifold-based mass lumping scheme.
   *
   * @param rowSumsHighOrder Row sums of the consistent mass matrix built from the element's own
   *        (high-order) shape functions, i.e. the integral of each \f$N_i\f$ over the element.
   * @param rowSumsLinear Row sums built from the linear (corner-node) shape functions of the same
   *        element, one entry per corner node.
   * @return The blend weight \f$w\f$ to use in \f$\hat{N} = w\,N + (1-w)\,N_\mathrm{lin}\f$.
   *
   * @details The scheme is manifold-based lumping:
   *
   *  - Yang, Y., Zheng, H. & Sivaselvan, M. V. (2017). "A rigorous and unified mass lumping scheme
   *    for higher-order elements". CMAME 319, 491-514. https://doi.org/10.1016/j.cma.2017.03.011
   *  - Duczek, S. & Gravenkamp, H. (2019). "Critical assessment of different mass lumping schemes
   *    for higher order serendipity finite elements". CMAME 350, 836-897.
   *    https://doi.org/10.1016/j.cma.2019.03.028
   *
   * It row-sums a blended shape function \f$\hat{N} = w N + (1-w) N_\mathrm{lin}\f$. Row-summing
   * \f$N\f$ alone is not positive for a serendipity element -- its corner shape functions integrate
   * to a negative value -- and blending in the linear shape functions is what restores positivity.
   *
   * The weight cannot be a constant, which is the subtlety this function exists for. A corner node's
   * lumped mass is \f$w S^{N}_i + (1-w) S^{\mathrm{lin}}_i\f$, so with \f$S^{N}_i < 0 < S^{\mathrm{lin}}_i\f$
   * positivity requires
   * \f[ w < w_\mathrm{max} = \min_i \frac{S^{\mathrm{lin}}_i}{S^{\mathrm{lin}}_i - S^{N}_i}. \f]
   * That limit is element-dependent, because the corner row sum is:
   *
   *  - quad8  (2D serendipity): \f$S^N = -\tfrac{1}{3}\f$ per corner (area 4)   -> \f$w_\mathrm{max} = 0.75\f$
   *  - hexa20 (3D serendipity): \f$S^N = -1\f$ per corner (volume 8)            -> \f$w_\mathrm{max} = 0.50\f$
   *
   * A hard-coded \f$w = \tfrac{1}{2}\f$ therefore sits comfortably inside the admissible range in 2D
   * and **exactly on its boundary** in 3D, where it produces corner masses of exactly zero -- a
   * singular lumped mass matrix, and one that lands a few 1e-17 either side of zero in floating point
   * rather than on it, so an "is any mass zero" check downstream passes and the inverse mass comes
   * out around 1e17.
   *
   * Returned here is \f$w = \min(\tfrac{1}{2}, \tfrac{2}{3} w_\mathrm{max})\f$. The
   * \f$\tfrac{2}{3}\f$ is a safety fraction chosen here rather than taken from either reference,
   * but it is not free-floating: it is the largest simple fraction that still yields
   * \f$w = \tfrac{1}{2}\f$ for quad8, and the values it produces are the published ones on both
   * element types where a published value exists --
   *
   *  - \f$w = \tfrac{1}{2}\f$ for quad8, the split of Yang et al. (2017),
   *  - \f$w = \tfrac{1}{3}\f$ for hexa20, the split assessed by Duczek & Gravenkamp (2019),
   *
   * so the rule is not merely a way of staying off the positivity boundary; it reproduces the
   * literature wherever the literature has an answer, and extrapolates in the same spirit where it
   * does not. The cap keeps \f$\tfrac{1}{2}\f$ as the value used wherever it is safe.
   *
   * An element with no negatively-integrating corner shape function (every linear element, since
   * there \f$N_\mathrm{lin} \equiv N\f$ and the result does not depend on \f$w\f$ at all) is
   * unconstrained and gets \f$\tfrac{1}{2}\f$.
   *
   * @note The reproduction of \f$\tfrac{1}{2}\f$ for quad8 is exact analytically, not necessarily
   * bit-exact: the row sums are accumulated numerically over the element's own quadrature points,
   * so \f$w_\mathrm{max}\f$ can land an ulp below \f$\tfrac{3}{4}\f$ (it does for a quad8 under
   * 2x2 reduced integration) and the returned weight then sits an ulp below \f$\tfrac{1}{2}\f$
   * rather than on it. The difference is far below any tolerance the lumped masses are checked
   * against, but it is a rounding-level difference and not a bitwise identity.
   *
   * @note Positivity is enforced at the corner nodes only, which is all a blend weight can reach.
   * A node with no linear counterpart carries \f$w S^{N}_i\f$, so if its own shape function
   * integrates to a negative value -- not possible for any of the element types in use, but
   * reachable in principle by a sufficiently distorted one -- its lumped mass is negative for
   * every \f$w > 0\f$ and this function cannot help. A positive returned weight is therefore not
   * a guarantee of a positive definite lumped mass matrix.
   */
  inline double manifoldBlendWeight( const Eigen::VectorXd& rowSumsHighOrder, const Eigen::VectorXd& rowSumsLinear )
  {
    constexpr double documentedWeight = 0.5;
    constexpr double fractionOfLimit  = 2.0 / 3.0;

    double admissibleWeight = std::numeric_limits< double >::max();
    for ( Eigen::Index i = 0; i < rowSumsLinear.size(); i++ ) {
      if ( rowSumsHighOrder( i ) >= 0.0 )
        continue;
      const double limit = rowSumsLinear( i ) / ( rowSumsLinear( i ) - rowSumsHighOrder( i ) );
      admissibleWeight   = std::min( admissibleWeight, limit );
    }

    if ( admissibleWeight == std::numeric_limits< double >::max() )
      return documentedWeight;

    return std::min( documentedWeight, fractionOfLimit * admissibleWeight );
  }

  /**
   * @brief The per-node lumped mass fractions of the manifold-based scheme, summing to one.
   *
   * @param rowSumsHighOrder Row sums from the element's own shape functions.
   * @param rowSumsLinear Row sums from the linear shape functions, one per corner node.
   * @param weight The blend weight, from manifoldBlendWeight().
   * @return One fraction per node, summing to one.
   *
   * @details Density-free on purpose: positivity and the mass *distribution* are properties of the
   * element geometry and its shape functions, so the same fractions serve the lumped mass (scaled by
   * density) and the critical time step (which needs only the distribution). Sharing them is what
   * keeps the two from drifting apart -- a time step derived from a different mass distribution than
   * the one actually assembled is exactly the sort of inconsistency that shows up as an unexplained
   * instability.
   *
   * Note the total is independent of the weight: the blend moves mass between the corner and the
   * remaining nodes but conserves the element total. That is why an incorrect weight leaves the
   * element mass, and hence the model mass, perfectly correct -- and is invisible to any check that
   * looks at totals.
   *
   * @note Integrating the blended shape function once, as done here and by the callers, is exactly
   * equivalent to assembling the blended consistent mass matrix
   * \f$M = \int \hat{N}^T \hat{N} \rho \,dV\f$ and row-summing that -- the form in which the scheme
   * is usually written and implemented. Both \f$N\f$ and \f$N_\mathrm{lin}\f$ are partitions of
   * unity and the blend is convex, so \f$\sum_j \hat{N}_j \equiv 1\f$ and
   * \f[ \sum_j M_{ij} = \int \hat{N}_i \Big( \sum_j \hat{N}_j \Big) \rho \,dV
   *                   = \int \hat{N}_i \rho \,dV. \f]
   * The linear form is therefore not an approximation of the quadratic one; it is the same numbers
   * at a fraction of the cost, and no consistent mass matrix has to be formed to get them.
   *
   * @note The normalisation is skipped, and the raw blended row sums returned, if their total is
   * not positive. That cannot happen for any element whose geometry is valid, since the total is
   * the element volume; it is a guard against dividing by a degenerate element's zero volume, not
   * a supported mode.
   */
  inline Eigen::VectorXd manifoldMassFractions( const Eigen::VectorXd& rowSumsHighOrder,
                                                const Eigen::VectorXd& rowSumsLinear,
                                                double                 weight )
  {
    Eigen::VectorXd fractions = weight * rowSumsHighOrder;
    fractions.head( rowSumsLinear.size() ) += ( 1.0 - weight ) * rowSumsLinear;

    const double total = fractions.sum();
    if ( total > 0.0 )
      fractions /= total;

    return fractions;
  }

  /**
   * @brief The factor by which a non-uniform lumped mass tightens the critical time step.
   *
   * @param massFractions Per-node lumped mass fractions, from manifoldMassFractions().
   * @return The factor to apply to a characteristic-length / wave-speed time step estimate.
   *
   * @details The usual estimate \f$\Delta t = l / c\f$ is the stable increment for an element whose
   * mass is spread uniformly over its nodes, each carrying \f$1/n\f$ of it. Lumping does not spread
   * it uniformly, and the lightest node sets the highest frequency: \f$\omega = \sqrt{k/m}\f$, so the
   * stable increment scales with \f$\sqrt{m_\mathrm{min} / m_\mathrm{uniform}}\f$, which is what is
   * returned.
   *
   * Exactly 1 for a linear element on a regular mesh, slightly below 1 for a distorted one (correct:
   * a distorted element does have a tighter limit than its volume alone suggests), and 0.913 for a
   * hexa20 under this scheme.
   *
   * @note This corrects only the mass-distribution part of the estimate. It does NOT correct for
   * polynomial order, and that residue is large: a convergence study on a GC3D20R bar shows the
   * answer settling only around a courant number of 0.1-0.2, i.e. roughly a further factor of five
   * unaccounted for, because \f$l/c\f$ is a linear-element formula and a quadratic element's highest
   * free eigenfrequency is well above what it predicts. Closing that properly wants an
   * eigenvalue-based estimate rather than another factor.
   *
   * @warning A non-positive smallest fraction returns 1.0, i.e. the uncorrected \f$l/c\f$
   * estimate. That is deliberately a fallback and not a fix: a non-positive nodal mass means the
   * lumping itself has failed, no time step derived from it is meaningful, and no blend weight can
   * repair it (see the second @note on manifoldBlendWeight()). Scaling \f$l/c\f$ by
   * \f$\sqrt{m_\mathrm{min}/m_\mathrm{uniform}}\f$ would return zero or NaN and mask the cause, so
   * the estimate is left alone and the assembled mass is where the failure belongs.
   *
   * @warning Nothing currently reports it there either: the only check on the assembled lumped
   * mass is on the solver side and tests for an exact zero, which a negative or near-zero nodal
   * mass passes. So a lumping failure is at present silent in both places. Detecting it at
   * assembly is the right fix and is deliberately out of scope here.
   */
  inline double timeStepFactorFromMassDistribution( const Eigen::VectorXd& massFractions )
  {
    if ( massFractions.size() == 0 )
      return 1.0;

    const double smallest = massFractions.minCoeff();
    if ( smallest <= 0.0 )
      return 1.0;

    return std::sqrt( static_cast< double >( massFractions.size() ) * smallest );
  }

} // namespace Marmot::FiniteElement::MassLumping
