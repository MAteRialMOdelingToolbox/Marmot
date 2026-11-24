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
#include "Marmot/MarmotFiniteElement.h"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Marmot::Discretizations::Elements::MassLumping {

  /**
   * @brief The number of nodes of the linear element whose shape functions are blended in.
   *
   * @tparam nDim   Number of spatial dimensions.
   * @tparam nNodes Number of element nodes.
   * @return The number of corner nodes of the reference domain of the element: e.g. 2 for bar3, 3 for tri3, 4 for
   *         quad8 and tetra10, 8 for hexa20. A linear element is its own linear counterpart.
   *
   * @details The blend adds the row sums of the linear element to the first nodes of the element, so the corner
   * nodes must come first in the node order, as in the Abaqus convention. An unknown combination of @p nDim and
   * @p nNodes is a compile error.
   */
  template < int nDim, int nNodes >
  constexpr int nLinearNodes()
  {
    return getNumberOfCornerNodes( getReferenceDomain( getElementShapeByMetric( nDim, nNodes ) ) );
  }

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

} // namespace Marmot::Discretizations::Elements::MassLumping
