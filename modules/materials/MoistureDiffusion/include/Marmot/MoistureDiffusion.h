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
#include "Marmot/MarmotTransportMaterial.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::MoistureDiffusion
   * @brief Moisture diffusion in concrete in terms of the pore relative humidity: Bažant–Najjar diffusivity, a
   * sorption isotherm for the moisture capacity, and an optional self-desiccation sink.
   *
   * One field, the pore relative humidity \f$h\f$ (0 ... 1). The content is the evaporable water content
   * \f$w(h)\f$ (sorption isotherm), the flux and the source are
   * \f[
   *   \mathbf{q} = -D(h)\,\nabla h,\qquad
   *   D(h) = D_1\left(\alpha_0 + \frac{1 - \alpha_0}{1 + \left(\frac{1-h}{1-h_c}\right)^n}\right),\qquad
   *   r(t) = -\dot w_s(t),
   * \f]
   * so that the balance reads \f$\frac{\partial w}{\partial h}\dot h - \nabla\cdot(D(h)\nabla h) = -\dot w_s\f$
   * (Bažant & Najjar, Matériaux et Constructions 5 (1972) 3–20). For \f$h \ge 1\f$ the diffusivity is held at
   * \f$D_1\f$.
   *
   * Sorption isotherms \f$w(h)\f$:
   * - type 0, linear: \f$w = K h\f$, the constant moisture capacity \f$K = \partial w/\partial h\f$. \f$K = 1\f$ gives
   *   the classical form \f$\dot h = \nabla\cdot(C(h)\nabla h)\f$, with \f$D_1 = C_1\f$ the diffusivity in m²/s.
   *   \f$K = 0\f$: no storage, steady state diffusion.
   * - type 1, BSB (Brunauer, Skalny & Bodor, J. Colloid Interface Sci. 30 (1969) 546–552; for cementitious materials
   *   Xi, Bažant & Jennings, Adv. Cem. Based Mater. 1 (1994) 248–257):
   *   \f$w = \dfrac{C\,k\,V_m\,h}{(1 - k h)(1 + (C - 1) k h)}\f$, with \f$C > 0\f$, \f$0 < k < 1\f$ and the monolayer
   *   capacity \f$V_m > 0\f$.
   *
   * Self-desiccation (optional): the water consumed by hydration in the sealed state,
   * \f$w_s(t) = \Delta w_\infty (1 - e^{-t/\tau_s})\f$, with the analysis time \f$t\f$ taken as the time since the
   * onset of self-desiccation. Over an increment \f$[t - \Delta t, t]\f$ the sink is the exact average rate,
   * \f$r = -\Delta w_\infty (e^{-(t - \Delta t)/\tau_s} - e^{-t/\tau_s})/\Delta t\f$, so that a sealed specimen
   * follows \f$w(h(t)) = w(h_0) - w_s(t)\f$ exactly at the time steps. It does not depend on \f$h\f$.
   *
   * Units: \f$w\f$, \f$K\f$, \f$V_m\f$, \f$\Delta w_\infty\f$ in a content unit, e.g. kg/m³; then \f$D_1\f$ in
   * kg/(m s) and a film coefficient in kg/(m² s). With \f$K = 1\f$, everything is in units of \f$h\f$: \f$D_1\f$ in
   * m²/s, film coefficients in m/s.
   *
   * Material properties, in this order (4, 6, 8 or 10 of them):
   *
   * | # | symbol | meaning | unit (SI) |
   * |---|--------|---------|-----------|
   * | 0 | \f$D_1\f$ | diffusivity at \f$h = 1\f$ | content m²/s |
   * | 1 | \f$\alpha_0\f$ | ratio \f$D(0)/D_1\f$ (approximately) | – |
   * | 2 | \f$h_c\f$ | relative humidity at the drop of \f$D\f$ (\f$D(h_c) = D_1 (1 + \alpha_0)/2\f$) | – |
   * | 3 | \f$n\f$ | steepness of the drop, \f$n \ge 1\f$ | – |
   * | 4 | type | sorption isotherm: 0 linear, 1 BSB | – |
   * | 5 | \f$K\f$ (type 0) or \f$C\f$ (type 1) | moisture capacity, or BSB energy constant | content, – |
   * | 6 | \f$k\f$ (type 1) | BSB constant | – |
   * | 7 | \f$V_m\f$ (type 1) | BSB monolayer capacity | content |
   * | 8 | \f$\Delta w_\infty\f$ | final self-desiccation water consumption, 0: none | content |
   * | 9 | \f$\tau_s\f$ | self-desiccation time constant | s |
   *
   * Short forms: 4 properties (linear isotherm, \f$K = 1\f$, no sink), 6 properties (type 0 with \f$K\f$), 8 properties
   * (no sink). Entries 6 and 7 are not used by the linear isotherm.
   *
   * @throws std::invalid_argument for another number of properties or parameters out of range.
   */
  class MoistureDiffusion : public TransportMaterial {
  public:
    /// The sorption isotherms
    enum IsothermType {
      Linear = 0,
      BSB    = 1,
    };

    /**
     * @brief Construct the material.
     * @param materialProperties Pointer to the material properties, see the class description.
     * @param nMaterialProperties Number of material properties: 4, 6, 8 or 10.
     * @param materialNumber Material identifier.
     */
    MoistureDiffusion( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /// @return 1, the relative humidity
    int getNumberOfFields() const override { return 1; }

    /// @return False for a linear isotherm with \f$K = 0\f$: steady state diffusion
    bool hasStorage( [[maybe_unused]] int field ) const override { return !( isotherm == Linear && K == 0.0 ); }

    /**
     * @brief Water content increment, moisture flux, self-desiccation sink and their algorithmic tangent.
     * @param[inout] response The response.
     * @param[out] tangent The algorithmic tangent.
     * @param[in] fieldState Relative humidities and gradient.
     * @param[in] timeIncrement The increment [time - dT, time] (self-desiccation).
     * @throws std::runtime_error if the BSB isotherm is evaluated outside its range (\f$k h \ge 1\f$ or
     * \f$1 + (C-1) k h \le 0\f$).
     */
    void computeTransport( Response&            response,
                           AlgorithmicTangent&  tangent,
                           const FieldState&    fieldState,
                           const TimeIncrement& timeIncrement ) const override;

    /**
     * @brief The Bažant–Najjar diffusivity.
     * @param h The relative humidity.
     * @param[out] dD_dh Its derivative.
     * @return \f$D(h)\f$.
     */
    double diffusivity( double h, double& dD_dh ) const;

    /**
     * @brief The sorption isotherm.
     * @param h The relative humidity.
     * @param[out] dw_dh The moisture capacity \f$\partial w/\partial h\f$.
     * @return \f$w(h)\f$.
     */
    double waterContent( double h, double& dw_dh ) const;

    /**
     * @brief The self-desiccation source, the average rate over the increment.
     * @param time The time at the end of the increment.
     * @param dT The time increment.
     * @return \f$r = -\dot w_s\f$ (non-positive).
     */
    double selfDesiccationSource( double time, double dT ) const;

  private:
    double       D1;     ///< \f$D_1\f$
    double       alpha0; ///< \f$\alpha_0\f$
    double       hc;     ///< \f$h_c\f$
    double       n;      ///< \f$n\f$
    IsothermType isotherm;
    double       K;      ///< linear isotherm: capacity
    double       C;      ///< BSB: C
    double       k;      ///< BSB: k
    double       Vm;     ///< BSB: monolayer capacity
    double       dwInf;  ///< self-desiccation: final water consumption
    double       tauS;   ///< self-desiccation: time constant
  };

} // namespace Marmot::Materials
