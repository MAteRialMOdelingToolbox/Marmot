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

namespace Marmot::ContinuumMechanics::TimeIntegration {

  /**
   * @brief Newmark-beta time integration at one point: the acceleration and the velocity at the end of an increment
   * from the displacement increment.
   * @details Solves the Newmark relations
   * @f[
   *   \Delta u = \Delta t\, v_n + \tfrac{\Delta t^2}{2} \left[ (1 - 2\beta)\, a_n + 2\beta\, a_{n+1} \right],
   *   \qquad
   *   v_{n+1} = v_n + \Delta t \left[ (1 - \gamma)\, a_n + \gamma\, a_{n+1} \right]
   * @f]
   * for @f$ a_{n+1} @f$ and @f$ v_{n+1} @f$, component by component, and returns the derivative
   * @f$ \partial a_{n+1} / \partial \Delta u = 1 / (\beta \Delta t^2) @f$, which the host's tangent needs for the
   * inertia term.
   *
   * Two cases add no inertia:
   * - @f$ \beta = 0 @f$: no time integration (a quasi-static analysis), which is at rest: the velocity and the
   *   acceleration are zero, whatever they were before (e.g. after a dynamic step), and so is the derivative.
   * - A zero-length increment, @f$ \Delta t = 0 @f$: nothing moves in zero time, so the velocity and the acceleration
   *   are left as they are (@f$ v_{n+1} = v_n @f$, @f$ a_{n+1} = a_n @f$), whatever @f$ \Delta u @f$ is, and the
   *   derivative is zero. The inertial force of the unchanged acceleration remains in the host's residual. A host
   *   that integrates in time skips such an increment, or evaluates at it only to compute a state at rest in time
   *   (e.g. the acceleration of its own nodes from equilibrium).
   *
   * @tparam nDim Number of components.
   * @param[in] du Displacement increment @f$ \Delta u @f$, nDim components.
   * @param[in,out] v Velocity: @f$ v_n @f$ on entry, @f$ v_{n+1} @f$ on return.
   * @param[in,out] a Acceleration: @f$ a_n @f$ on entry, @f$ a_{n+1} @f$ on return.
   * @param[in] dT Time increment @f$ \Delta t \geq 0 @f$.
   * @param[in] newmarkBeta Newmark parameter @f$ \beta @f$ (1/4 with @f$ \gamma = 1/2 @f$: the average acceleration
   * rule; 0: no time integration).
   * @param[in] newmarkGamma Newmark parameter @f$ \gamma @f$.
   * @param[out] da_ddu Derivative @f$ \partial a_{n+1} / \partial \Delta u @f$, an nDim x nDim matrix (column
   * major); only its diagonal is written.
   */
  template < int nDim >
  void newmarkBetaIntegration( const double* du,
                               double*       v,
                               double*       a,
                               double        dT,
                               double        newmarkBeta,
                               double        newmarkGamma,
                               double*       da_ddu )
  {
    if ( newmarkBeta == 0.0 ) {
      for ( int i = 0; i < nDim; i++ ) {
        v[i]                 = 0.0;
        a[i]                 = 0.0;
        da_ddu[i + i * nDim] = 0.0;
      }
      return;
    }

    if ( dT == 0.0 ) {
      for ( int i = 0; i < nDim; i++ )
        da_ddu[i + i * nDim] = 0.0;
      return;
    }

    for ( int i = 0; i < nDim; i++ ) {
      const double du_tilde = dT * v[i] + 0.5 * dT * dT * ( ( 1 - 2 * newmarkBeta ) * a[i] );
      const double v_tilde  = v[i] + dT * ( 1 - newmarkGamma ) * a[i];

      a[i]                 = ( du[i] - du_tilde ) / ( newmarkBeta * dT * dT );
      da_ddu[i + i * nDim] = 1. / ( newmarkBeta * dT * dT );
      v[i]                 = v_tilde + newmarkGamma * dT * a[i];
    }
  }

} // namespace Marmot::ContinuumMechanics::TimeIntegration
