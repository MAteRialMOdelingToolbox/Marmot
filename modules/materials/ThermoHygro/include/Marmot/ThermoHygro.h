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
#include "Marmot/HeatConduction.h"
#include "Marmot/MarmotTransportMaterial.h"
#include "Marmot/MoistureDiffusion.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::ThermoHygro
   * @brief Coupled heat and moisture transport in concrete: the temperature \f$T\f$ and the pore relative humidity
   * \f$h\f$, Bažant–Thonguthai type.
   *
   * Two fields, in this order: 0 the temperature \f$T\f$, 1 the pore relative humidity \f$h\f$. The material extends
   * Marmot::Materials::HeatConduction (properties 0–10) and Marmot::Materials::MoistureDiffusion (properties 11–20),
   * which it uses for the uncoupled parts, by the couplings (properties 21–31). With all couplings off, both fields
   * follow the two uncoupled materials exactly. With \f$\theta = T - T_\mathrm{ref}\f$, the absolute temperature
   * \f$T_a = T - T_0\f$, and \f$\hat h = h - \kappa\,\theta\f$:
   *
   * - water content (temperature dependent sorption isotherm, hygrothermic coefficient \f$\kappa =
   *   (\partial h/\partial T)_w\f$): \f$w(h, T) = w_\mathrm{iso}(\hat h)\f$, \f$w_\mathrm{iso}\f$ the isotherm of
   *   MoistureDiffusion;
   * - moisture flux, driven by \f$\nabla h\f$ and, for \f$\chi = 1\f$, by the gradient of the vapour pressure
   *   \f$p_v = h\,p_s(T)\f$ (Bažant–Thonguthai): \f$\mathbf{q}_w = -D(h,T)\,(\nabla h + \chi\,h\,g(T)\,\nabla T)\f$,
   *   \f$g = \mathrm{d}\ln p_s/\mathrm{d}T = \Theta_v/T_a^2\f$ (Clausius–Clapeyron, \f$\Theta_v = M_w L_v/R\f$),
   *   \f$D(h,T) = D_\mathrm{BN}(h)\,f(T)\f$, \f$f(T) = (T_a/T_{a,\mathrm{ref}})^{m}
   *   \exp(Q_D/R\,(1/T_{a,\mathrm{ref}} - 1/T_a))\f$ (Arrhenius);
   * - enthalpy increment, with the sensible heat of the evaporable water:
   *   \f$\Delta H = H_\mathrm{dry}(T) - H_\mathrm{dry}(T_n) + c_w\,w(h,T)\,(T - T_n)\f$;
   * - heat flux, moisture dependent conductivity and latent heat carried by the vapour (Künzel):
   *   \f$\mathbf{q}_T = -(1 + \beta_{kh}h)\,\mathbf{k}(T)\,\nabla T + \ell\,\mathbf{q}_w\f$, \f$\ell = \psi L_v\f$ with
   *   the vapour fraction \f$\psi\f$ of the moisture flux;
   * - sources, as functions of the equivalent (maturity) age \f$t_e\f$, \f$\dot t_e = \exp(E_h/R\,(1/T_{a,\mathrm{ref}}
   *   - 1/T_a))\f$: heat of hydration \f$Q(t_e) = Q_\infty(1 - e^{-t_e/\tau_q})\f$, \f$r_T = \dot Q\f$, and
   *   self-desiccation \f$w_s(t_e) = \Delta w_\infty(1 - e^{-t_e/\tau_s})\f$, \f$r_h = -\dot w_s\f$; both are averaged
   *   exactly over the increment, with \f$\dot t_e\f$ at the end of the increment (backward Euler). For \f$E_h = 0\f$,
   *   \f$t_e\f$ is the analysis time, as in MoistureDiffusion.
   *
   * Material properties, in this order (21 or 32 of them):
   *
   * | # | symbol | meaning | unit (SI) |
   * |---|--------|---------|-----------|
   * | 0–10 | | HeatConduction: \f$\rho, c_0, k_{11}, k_{22}, k_{33}, k_{12}, k_{13}, k_{23}, \beta_k, \beta_c,
   *   T_\mathrm{ref}\f$ | |
   * | 11–20 | | MoistureDiffusion: \f$D_1, \alpha_0, h_c, n\f$, isotherm type, \f$K\f$ or \f$C\f$, \f$k\f$, \f$V_m\f$,
   *   \f$\Delta w_\infty\f$, \f$\tau_s\f$ | |
   * | 21 | \f$T_0\f$ | absolute zero of the temperature scale: -273.15 (°C) or 0 (K) | K or °C |
   * | 22 | \f$Q_D/R\f$ | activation temperature of the diffusivity, 0: none | K |
   * | 23 | \f$m\f$ | exponent of \f$T_a/T_{a,\mathrm{ref}}\f$ in \f$f(T)\f$ (Bažant–Najjar: 1) | – |
   * | 24 | \f$\kappa\f$ | hygrothermic coefficient \f$(\partial h/\partial T)_w\f$ | 1/K |
   * | 25 | \f$\chi\f$ | 0: flux driven by \f$\nabla h\f$ only, 1: by the vapour pressure gradient | – |
   * | 26 | \f$c_w\f$ | specific heat of the evaporable water | J/(content K) |
   * | 27 | \f$\beta_{kh}\f$ | moisture coefficient of the conductivity | – |
   * | 28 | \f$\ell\f$ | latent heat carried per unit moisture flux, \f$\psi L_v\f$ | J/content |
   * | 29 | \f$E_h/R\f$ | activation temperature of the hydration (equivalent age), 0: none | K |
   * | 30 | \f$Q_\infty\f$ | final heat of hydration, 0: none | J/m³ |
   * | 31 | \f$\tau_q\f$ | time constant of the heat of hydration | s |
   *
   * 21 properties: no coupling. The full forms of both uncoupled materials are required (11 and 10 properties).
   * The content unit is that of the isotherm (e.g. kg/m³).
   *
   * State variables: "equivalent age" \f$t_e\f$ (at the end of the last converged increment; starts at 0).
   *
   * @throws std::invalid_argument for another number of properties or parameters out of range.
   */
  class ThermoHygro : public TransportMaterial {
  public:
    /// \f$\Theta_v = M_w L_v / R\f$ = 0.01801528 kg/mol * 2.45e6 J/kg / 8.314462618 J/(mol K), Clausius–Clapeyron
    static constexpr double vapourTemperature = 0.01801528 * 2.45e6 / 8.314462618;

    /**
     * @brief Construct the material.
     * @param materialProperties Pointer to the material properties, see the class description.
     * @param nMaterialProperties Number of material properties: 21 or 32.
     * @param materialNumber Material identifier.
     */
    ThermoHygro( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /// @return 2: the temperature and the relative humidity
    int getNumberOfFields() const override { return 2; }

    /// @return Field 0 (temperature): heat capacity or water with heat capacity; field 1: as MoistureDiffusion
    bool hasStorage( int field ) const override;

    /**
     * @brief Content increments, fluxes and sources of heat and moisture, and their algorithmic tangent (all four
     * blocks).
     * @param[inout] response The response; state variable "equivalent age" updated.
     * @param[out] tangent The algorithmic tangent.
     * @param[in] fieldState Temperatures, relative humidities and their gradients.
     * @param[in] timeIncrement The increment [time - dT, time].
     * @throws std::runtime_error outside the range of the laws (non-positive absolute temperature or conductivity
     * factor, BSB isotherm out of range).
     */
    void computeTransport( Response&            response,
                           AlgorithmicTangent&  tangent,
                           const FieldState&    fieldState,
                           const TimeIncrement& timeIncrement ) const override;

    /**
     * @brief The temperature factor of the diffusivity.
     * @param T The temperature.
     * @param[out] df_dT Its derivative.
     * @return \f$f(T)\f$.
     */
    double diffusivityTemperatureFactor( double T, double& df_dT ) const;

    /**
     * @brief The rate of the equivalent age (maturity), \f$\exp(E_h/R\,(1/T_{a,\mathrm{ref}} - 1/T_a))\f$.
     * @param T The temperature.
     * @param[out] da_dT Its derivative.
     * @return The rate.
     */
    double maturityRate( double T, double& da_dT ) const;

  private:
    HeatConduction    heat;     ///< the uncoupled heat conduction law
    MoistureDiffusion moisture; ///< the uncoupled moisture diffusion law

    double TReference;          ///< \f$T_\mathrm{ref}\f$
    double dwInf;               ///< \f$\Delta w_\infty\f$
    double tauS;                ///< \f$\tau_s\f$
    double TZero;               ///< \f$T_0\f$
    double QD;                  ///< \f$Q_D/R\f$
    double mD;                  ///< \f$m\f$
    double kappa;               ///< \f$\kappa\f$
    double chi;                 ///< \f$\chi\f$
    double cw;                  ///< \f$c_w\f$
    double betaKh;              ///< \f$\beta_{kh}\f$
    double ell;                 ///< \f$\ell\f$
    double Eh;                  ///< \f$E_h/R\f$
    double QInf;                ///< \f$Q_\infty\f$
    double tauQ;                ///< \f$\tau_q\f$
    bool   withMaturity;        ///< \f$E_h \neq 0\f$: equivalent age state

    /// absolute temperature, throws if not positive
    double absoluteTemperature( double T ) const;
  };

} // namespace Marmot::Materials
