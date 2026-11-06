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
#include "Marmot/MarmotBondSlipLaw.h"
#include <utility>

namespace Marmot::Materials {

  /**
   * @class ModelCode2010BondSlip
   * @brief Bond-slip law with the monotonic envelope of the fib Model Code 2010 (Sec. 6.1.1) and elastic unloading.
   *
   * Monotonic envelope of the tangential bond stress \f$\tau_t\f$ over the slip \f$s \geq 0\f$:
   * \f[
   *   \tau_{\mathrm{env}}(s) = \begin{cases}
   *     K_0\, s & 0 \leq s \leq s_0 \\
   *     \tau_{\max} \left( s / s_1 \right)^{\alpha} & s_0 < s \leq s_1 \\
   *     \tau_{\max} & s_1 < s \leq s_2 \\
   *     \tau_{\max} - (\tau_{\max} - \tau_f) \frac{s - s_2}{s_3 - s_2} & s_2 < s \leq s_3 \\
   *     \tau_f & s_3 < s
   *   \end{cases}
   * \f]
   * The linear branch with the stiffness \f$K_0\f$ regularizes the infinite initial slope of the power law; \f$s_0\f$
   * is where both meet.
   *
   * Unloading and reloading are elastic with \f$K_0\f$, \f$\tau_t = K_0 (s - s_p)\f$, bounded by
   * \f$|\tau_t| \leq \tau_{\mathrm{env}}(\kappa)\f$, where \f$\kappa = \max |s|\f$ is the largest slip ever reached
   * in either direction. Beyond that bound the slip is frictional (plastic slip \f$s_p\f$). The bond damage is
   * therefore symmetric, the reversed direction is bounded by the same, already reduced bond stress. Cyclic
   * degradation beyond that is not modelled.
   *
   * The normal behavior is linear elastic with the stiffness \f$K_n\f$.
   *
   * Material properties:
   *  0. @b tauMax: bond strength \f$\tau_{\max}\f$
   *  1. @b s1, 2. @b s2, 3. @b s3: characteristic slips, \f$0 < s_1 \leq s_2 < s_3\f$
   *  4. @b alpha: exponent of the ascending branch, \f$0 \leq \alpha < 1\f$
   *  5. @b tauF: residual (frictional) bond stress, \f$0 \leq \tau_f \leq \tau_{\max}\f$
   *  6. @b K0: initial and unloading bond stiffness, \f$K_0 > \tau_{\max} / s_1\f$
   *  7. @b Kn: normal stiffness
   *
   * State variables: @b plastic slip \f$s_p\f$, @b max slip \f$\kappa\f$.
   */
  class ModelCode2010BondSlip : public Marmot::Materials::BondSlipLaw {
  public:
    const double& tauMax; ///< bond strength
    const double& s1;     ///< slip at the begin of the plateau
    const double& s2;     ///< slip at the end of the plateau
    const double& s3;     ///< slip at the begin of the residual branch
    const double& alpha;  ///< exponent of the ascending branch
    const double& tauF;   ///< residual bond stress
    const double& K0;     ///< initial and unloading stiffness
    const double& Kn;     ///< normal stiffness
    const double  s0;     ///< end of the linear regularization of the ascending branch

    ModelCode2010BondSlip( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeBondStress( State&                        state,
                            Marmot::Math::Matrix3d&       dBondStress_dSlip,
                            const Marmot::Math::Vector3d& slip,
                            const Marmot::Math::Vector3d& dSlip,
                            const TimeIncrement&          timeIncrement ) const override;

    /**
     * @brief The monotonic envelope and its slope.
     * @param s slip, \f$s \geq 0\f$
     * @return \f$\{\tau_{\mathrm{env}}(s),\; \tau_{\mathrm{env}}'(s)\}\f$
     */
    std::pair< double, double > envelope( double s ) const;
  };
} // namespace Marmot::Materials
