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
#include <cstddef>
#include <limits>
#include <vector>

/**
 * @brief The quadrature of the consistent mass matrix, independent of the element's own rule.
 *
 * @details An element evaluates its stiffness and internal forces at the quadrature points of its
 * own rule, which for a reduced-integration element is one order below full. The consistent mass
 * \f$\mathbf{M}_e = \int \rho\, \mathbf{N}^\mathsf{T}\mathbf{N}\, dV\f$ cannot be integrated that
 * way: it has rank at most (number of quadrature points) \f$\times\f$ (number of components), so a
 * reduced 20-node hexahedron, with 8 points, yields a mass of rank 24 out of 60. In the dynamic
 * tangent \f$\mathbf{K} + \mathbf{M} / (\beta \Delta t^2)\f$ the stiffness hides that; wherever the
 * mass is solved on its own -- the initial acceleration \f$\mathbf{M}\,\mathbf{a}_0 = \mathbf{R}\f$
 * of an implicit dynamic analysis -- the system is singular, and the modes of the null space carry
 * no inertia at all.
 *
 * The consistent mass is therefore always integrated with the FULL Gauss rule of the element's
 * shape, whatever rule the element uses otherwise. For a full-integration element that is its own
 * rule, and the result is unchanged. The density (and any other inertia coefficient) is a material
 * quantity that lives at the element's own quadrature points; at a point of the mass rule it is
 * taken from the nearest of them in parent coordinates -- exact for a density that is constant
 * over the element, a piecewise-constant approximation otherwise.
 *
 * @note For the tetrahedra the "full" rules are the element rules themselves (1 point for Tetra4,
 * 4 for Tetra10), which are too low for a mass matrix as well. No tetrahedral element is currently
 * registered; one that is must be given a rule of sufficient order here first.
 */
namespace Marmot::FiniteElement::ConsistentMass {

  /**
   * @brief The quadrature rule the consistent mass of an element of the given shape is integrated
   * with: the full Gauss rule of that shape.
   * @param shape The element shape.
   * @return The quadrature points and weights.
   */
  inline const std::vector< Quadrature::QuadraturePointInfo >& integrationRule( ElementShapes shape )
  {
    return Quadrature::getGaussPointInfo( shape, Quadrature::IntegrationTypes::FullIntegration );
  }

  /**
   * @brief The element quadrature point nearest to a point of the mass rule, in parent coordinates.
   * @tparam QuadraturePointContainer A container of the element's quadrature points, each with a
   * member @c xi holding its parent coordinates.
   * @param xi Parent coordinates of the point of the mass rule.
   * @param qps The element's quadrature points.
   * @return The index of the nearest one; the first of them on a tie.
   */
  template < class QuadraturePointContainer >
  std::size_t nearestQuadraturePoint( const Eigen::VectorXd& xi, const QuadraturePointContainer& qps )
  {
    std::size_t nearest         = 0;
    double      nearestDistance = std::numeric_limits< double >::infinity();
    for ( std::size_t i = 0; i < qps.size(); i++ ) {
      const double distance = ( qps[i].xi - xi ).squaredNorm();
      if ( distance < nearestDistance ) {
        nearest         = i;
        nearestDistance = distance;
      }
    }
    return nearest;
  }

} // namespace Marmot::FiniteElement::ConsistentMass
