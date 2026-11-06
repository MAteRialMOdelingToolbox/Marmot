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
   * @class ModelCode2010CoulombBondSlip
   * @brief Bond-slip law with the fib Model Code 2010 envelope plus pressure-dependent (Coulomb) friction and a
   * limited normal stress (local crushing of the concrete).
   *
   * Normal behavior (components 1 and 2 of the local frame), elastic-perfectly plastic with an isotropic cap in the
   * normal plane: \f$\boldsymbol{t}_n = K_n (\boldsymbol{s}_n - \boldsymbol{s}_n^p)\f$,
   * \f$|\boldsymbol{t}_n| \leq p_{\max}\f$ (radial return; \f$p_{\max} \leq 0\f$: no cap). The plastic normal slip
   * \f$\boldsymbol{s}_n^p\f$ is the permanently crushed channel wall. The contact pressure is
   * \f$p = |\boldsymbol{t}_n|\f$: a bar pressed against its channel in any normal direction.
   *
   * Tangential behavior as ModelCode2010BondSlip (elastic with \f$K_0\f$, friction on a bound, history
   * \f$\kappa = \max |s|\f$), with the bound raised by Coulomb friction:
   * \f[ |\tau_t| \leq \tau_{\mathrm{env}}(\kappa) + \mu\, p . \f]
   * The tangent is consistent, including \f$\partial \tau_t / \partial \boldsymbol{s}_n = \mathrm{sign}(\tau_t)\,\mu\,
   * \partial p / \partial \boldsymbol{s}_n\f$ on the friction bound. For \f$\mu = 0\f$ and no cap, the law equals
   * ModelCode2010BondSlip.
   *
   * Material properties:
   *  0. @b tauMax, 1. @b s1, 2. @b s2, 3. @b s3, 4. @b alpha, 5. @b tauF, 6. @b K0, 7. @b Kn: as ModelCode2010BondSlip
   *  8. @b mu: friction coefficient, \f$\mu \geq 0\f$
   *  9. @b pMax: largest normal stress (crushing), \f$\leq 0\f$ for none
   *
   * State variables: @b plastic slip \f$s_p\f$, @b max slip \f$\kappa\f$, @b plastic normal slip (2).
   */
  class ModelCode2010CoulombBondSlip : public Marmot::Materials::BondSlipLaw {
  public:
    const double& tauMax; ///< bond strength
    const double& s1;     ///< slip at the begin of the plateau
    const double& s2;     ///< slip at the end of the plateau
    const double& s3;     ///< slip at the begin of the residual branch
    const double& alpha;  ///< exponent of the ascending branch
    const double& tauF;   ///< residual bond stress
    const double& K0;     ///< initial and unloading stiffness
    const double& Kn;     ///< normal stiffness
    const double& mu;     ///< friction coefficient
    const double& pMax;   ///< largest normal stress, <= 0: unlimited
    const double  s0;     ///< end of the linear regularization of the ascending branch

    ModelCode2010CoulombBondSlip( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeBondStress( State&                        state,
                            Marmot::Math::Matrix3d&       dBondStress_dSlip,
                            const Marmot::Math::Vector3d& slip,
                            const Marmot::Math::Vector3d& dSlip,
                            const TimeIncrement&          timeIncrement ) const override;

    /**
     * @brief The monotonic envelope (without friction) and its slope.
     * @param s slip, \f$s \geq 0\f$
     * @return \f$\{\tau_{\mathrm{env}}(s),\; \tau_{\mathrm{env}}'(s)\}\f$
     */
    std::pair< double, double > envelope( double s ) const;
  };
} // namespace Marmot::Materials
