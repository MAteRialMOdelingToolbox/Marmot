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
   * @brief Bond-slip law with the monotonic envelope of the fib Model Code 2010 (Sec. 6.1.1), elastic unloading
   * and optional pressure-dependent (Coulomb) friction with a limited normal stress.
   *
   * <b>From fib MC2010</b> (Eqs. 6.1-1 to 6.1-4) is only the monotonic envelope of the tangential bond stress
   * \f$\tau_t\f$ over the slip \f$s \geq 0\f$:
   * \f[
   *   \tau_{\mathrm{env}}(s) = \begin{cases}
   *     K_0\, s & 0 \leq s \leq s_0 \\
   *     \tau_{\max} \left( s / s_1 \right)^{\alpha} & s_0 < s \leq s_1 \\
   *     \tau_{\max} & s_1 < s \leq s_2 \\
   *     \tau_{\max} - (\tau_{\max} - \tau_f) \frac{s - s_2}{s_3 - s_2} & s_2 < s \leq s_3 \\
   *     \tau_f & s_3 < s
   *   \end{cases}
   * \f]
   * Parameters for the envelope are given in MC2010 Table 6.1-1, see the documentation page on bond-slip laws.
   *
   * <b>Modelling choices, not from MC2010:</b>
   *  - The linear branch with the stiffness \f$K_0\f$ regularizes the infinite initial slope of the power law;
   *    \f$s_0\f$ is where both meet (\f$K_0 s_0 = \tau_{\max} (s_0/s_1)^\alpha\f$; \f$s_0 = 0\f$ for \f$\alpha = 1\f$).
   *  - Unloading and reloading are elastic with \f$K_0\f$, \f$\tau_t = K_0 (s - s_p)\f$, inside the yield condition
   *    \f$|\tau_t| \leq \tau_{\mathrm{env}}(\kappa) + \mu p\f$, where \f$\kappa = \max |s|\f$ is the largest total
   *    slip ever reached in either direction. The bond damage is therefore symmetric: the reversed direction is
   *    bounded by the same, already reduced bond stress. Cyclic degradation beyond that is not modelled.
   *  - The normal behavior is elastic with the penalty stiffness \f$K_n\f$,
   *    \f$\boldsymbol{t}_n = K_n (\boldsymbol{s}_n - \boldsymbol{s}_n^p)\f$, optionally perfectly plastic with an
   *    isotropic cap \f$|\boldsymbol{t}_n| \leq p_{\max}\f$ in the normal plane (radial return). The cap models local
   *    crushing of the concrete around the bar: the plastic normal slip \f$\boldsymbol{s}_n^p\f$ is the permanent
   *    widening of the hole the bar sits in.
   *  - Coulomb friction \f$\mu p\f$ with the contact pressure \f$p = |\boldsymbol{t}_n|\f$ raises the bond bound.
   *    The magnitude is used because a bar pushed sideways presses against the concrete on one side of its hole,
   *    whatever the sign of the normal slip. MC2010 accounts for transverse pressure by its own factor
   *    \f$\Omega_p\f$ instead, which is not implemented.
   *
   * With 8 properties, \f$\mu = 0\f$ and the normal behavior is linear elastic.
   *
   * Material properties:
   *  0. @b tauMax: bond strength \f$\tau_{\max}\f$
   *  1. @b s1, 2. @b s2, 3. @b s3: characteristic slips, \f$0 < s_1 \leq s_2 < s_3\f$
   *  4. @b alpha: exponent of the ascending branch, \f$0 \leq \alpha \leq 1\f$
   *  5. @b tauF: residual bond stress, \f$0 \leq \tau_f \leq \tau_{\max}\f$
   *  6. @b K0: initial and unloading bond stiffness (penalty), \f$K_0 > \tau_{\max} / s_1\f$
   *  7. @b Kn: normal stiffness (penalty), \f$K_n \geq 0\f$
   *  8. @b mu (optional, default 0): friction coefficient, \f$\mu \geq 0\f$
   *  9. @b pMax (optional, default 0): largest normal stress, \f$\leq 0\f$ for none; \f$p_{\max} > 0\f$ requires
   *     \f$K_n > 0\f$
   *
   * State variables: @b plastic slip \f$s_p\f$, @b max slip \f$\kappa\f$, @b plastic normal slip
   * \f$\boldsymbol{s}_n^p\f$ (2).
   */
  class ModelCode2010BondSlip : public Marmot::Materials::BondSlipLaw {
  public:
    const double tauMax; ///< bond strength
    const double s1;     ///< slip at the begin of the plateau
    const double s2;     ///< slip at the end of the plateau
    const double s3;     ///< slip at the begin of the residual branch
    const double alpha;  ///< exponent of the ascending branch
    const double tauF;   ///< residual bond stress
    const double K0;     ///< initial and unloading stiffness
    const double Kn;     ///< normal stiffness
    const double mu;     ///< friction coefficient
    const double pMax;   ///< largest normal stress, <= 0: unlimited
    const double s0;     ///< end of the linear regularization of the ascending branch

    /**
     * @brief Construct the law and validate its properties.
     * @param materialProperties  the 8 to 10 properties listed in the class documentation
     * @param nMaterialProperties number of properties
     * @param materialNumber      identifier of this material
     * @throws std::invalid_argument for fewer than 8 or more than 10 properties, or invalid values
     */
    ModelCode2010BondSlip( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeBondStress( State&                        state,
                            Marmot::Math::Matrix3d&       dBondStress_dSlip,
                            const Marmot::Math::Vector3d& slip,
                            const Marmot::Math::Vector3d& dSlip,
                            const TimeIncrement&          timeIncrement ) const override;

    /**
     * @brief The monotonic envelope (without friction) and its slope.
     * @param s slip, \f$s \geq 0\f$
     * @return \f$\{\tau_{\mathrm{env}}(s),\; \mathrm{d}\tau_{\mathrm{env}} / \mathrm{d}s\}\f$
     */
    std::pair< double, double > envelope( double s ) const;
  };
} // namespace Marmot::Materials
