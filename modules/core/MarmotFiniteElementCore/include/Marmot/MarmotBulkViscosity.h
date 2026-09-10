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
#include <algorithm>
#include <cmath>

namespace Marmot::FiniteElement::BulkViscosity {

  /**
   * @brief Coefficients of the artificial bulk viscosity.
   *
   * @details Both are dimensionless. The defaults are the ones Abaqus/Explicit applies unless told
   * otherwise, and they are defaults rather than recommendations: \f$b_1\f$ is sized to damp the
   * highest resolvable frequency of the mesh, not to model any physical dissipation.
   */
  struct Coefficients {
    /** Linear coefficient \f$b_1\f$, active in compression and in expansion. */
    double linear = 0.0;
    /** Quadratic coefficient \f$b_2\f$, active in compression only. */
    double quadratic = 0.0;
    /**
     * Exponent \f$n\f$ of the optional degradation with the material's loss of stiffness; see
     * degradationFactor(). Zero -- the default -- disables it entirely, which is both the
     * pre-existing behaviour and the cheap one.
     */
    double degradation = 0.0;

    /** @brief Whether any bulk viscosity is active at all. */
    bool areActive() const { return linear != 0.0 || quadratic != 0.0; }

    /** @brief Whether the viscous stress is to be degraded with the material's stiffness. */
    bool isDegraded() const { return degradation > 0.0; }
  };

  /**
   * @brief Factor by which the viscous stress is degraded as the material loses stiffness.
   *
   * @param currentWaveSpeed Wave speed of the material at its CURRENT state, i.e. from the
   *        algorithmic tangent as it is now.
   * @param referenceWaveSpeed Wave speed of the same material in its virgin state.
   * @param exponent The exponent \f$n\f$. Zero returns exactly one.
   * @return \f$\left(c/c_0\right)^n\f$, clamped to \f$[0,1]\f$.
   *
   * @details Motivation. A fully cracked element carries no real stress, but the linear term of
   * the artificial viscosity does not know that and keeps transmitting
   * \f$b_1\rho c_0 L_e \dot{\varepsilon}_\mathrm{vol}\f$ across the crack. Since that term is
   * active in expansion, it acts as a viscous resistance to the crack OPENING which never relaxes,
   * and its work is charged to the fracture energy. Measured on a gradient-damage bar, artificial
   * bulk viscosity at the Abaqus/Explicit default coefficients inflated the dissipated energy per
   * unit fracture area by 27 to 55 percent depending on the mesh -- and applying it everywhere
   * EXCEPT the elements that damaged changed that energy by -1.6 percent. Essentially the whole
   * error is generated inside the damaged elements.
   *
   * This is the same argument that restricts the quadratic term to compression, applied to the
   * STATE instead of to the sign of the rate: the quadratic term is gated so that it cannot resist
   * a crack opening, and this gates the linear term so that it stops resisting one once the
   * material there has failed.
   *
   * What it actually measures. The ratio of wave speeds is the square root of the ratio of tangent
   * stiffnesses, so for a model whose tangent degrades as \f$(1-\omega)\,\mathbb{C}_0\f$:
   *
   *  - \f$n = 0\f$: off. No degradation, and the current wave speed is never evaluated.
   *  - \f$n = 1\f$: scales with the wave speed, i.e. with \f$\sqrt{1-\omega}\f$.
   *  - \f$n = 2\f$: scales with the tangent stiffness, i.e. with \f$1-\omega\f$ itself.
   *
   * @warning It degrades with the current TANGENT, not with a damage variable -- no material
   * interface here reports damage, and adding one would change every element's vtable. For a
   * quasi-brittle material in tension the two coincide, because the softening IS the damage. For a
   * model that merely yields, the algorithmic tangent also drops and the viscous stress is then
   * degraded by plastic flow rather than by cracking, which is not what this is for.
   *
   * @note Not free. The current wave speed comes from the material's algorithmic tangent, and
   * obtaining that costs a full constitutive evaluation per quadrature point per increment -- which
   * is precisely why the reference value is cached. That is the reason this is opt-in through a
   * named element property rather than the default.
   */
  inline double degradationFactor( double currentWaveSpeed, double referenceWaveSpeed, double exponent )
  {
    if ( exponent <= 0.0 || referenceWaveSpeed <= 0.0 )
      return 1.0;

    const double ratio = std::clamp( currentWaveSpeed / referenceWaveSpeed, 0.0, 1.0 );

    if ( exponent == 1.0 )
      return ratio;
    if ( exponent == 2.0 )
      return ratio * ratio;

    return std::pow( ratio, exponent );
  }

  /**
   * @brief The artificial viscous stress added to every normal stress component.
   *
   * @param volumetricStrainRate Trace of the strain rate, \f$\dot{\varepsilon}_\mathrm{vol} =
   *        \mathrm{tr}(\dot{\boldsymbol{\varepsilon}})\f$. Negative in compression.
   * @param density Current mass density \f$\rho\f$. If the run is mass scaled, this is the SCALED
   *        density, which is correct: the term exists to damp the modes of the system that is
   *        actually being integrated, and mass scaling changes those.
   * @param waveSpeed Dilatational wave speed \f$c_d\f$ of the material.
   * @param characteristicElementLength The element's smallest physical extent \f$L_e\f$.
   * @param coefficients The two dimensionless coefficients.
   * @return The scalar \f$\sigma_\mathrm{bv}\f$ to ADD to \f$\sigma_{xx}\f$, \f$\sigma_{yy}\f$ and
   *         \f$\sigma_{zz}\f$ (i.e. \f$\boldsymbol{\sigma}_\mathrm{bv} = \sigma_\mathrm{bv}\,
   *         \mathbf{I}\f$). Deviatoric components are untouched.
   *
   * @details
   * \f[
   *   \sigma_\mathrm{bv} = b_1\,\rho\,c_d\,L_e\,\dot{\varepsilon}_\mathrm{vol}
   *                      - \rho\,(b_2 L_e)^2\,\dot{\varepsilon}_\mathrm{vol}^2
   *                        \,H(-\dot{\varepsilon}_\mathrm{vol})
   * \f]
   *
   * following the standard element-level artificial viscosity of
   *
   *  - VonNeumann, J. & Richtmyer, R. D. (1950). "A method for the numerical calculation of
   *    hydrodynamic shocks". Journal of Applied Physics 21(3), 232-237.
   *    https://doi.org/10.1063/1.1699639
   *  - Landshoff, R. (1955). "A numerical method for treating fluid flow in the presence of
   *    shocks". Los Alamos Scientific Laboratory report LA-1930.
   *
   * as it is applied in explicit finite element codes (Abaqus/Explicit, LS-DYNA), where the linear
   * (Landshoff) term damps element-level ringing and the quadratic (Von Neumann-Richtmyer) term
   * spreads a shock over a few elements.
   *
   * **This is a numerical device, not a material model.** It adds no state, is not seen by the
   * constitutive law, and does not enter the stored stress: it is added to the stress only where
   * the internal force is integrated. Both terms are dissipative by construction -- the power
   * density is
   * \f$\sigma_\mathrm{bv}\dot{\varepsilon}_\mathrm{vol}
   *  = b_1\rho c_d L_e \dot{\varepsilon}_\mathrm{vol}^2
   *  - \rho (b_2 L_e)^2 \dot{\varepsilon}_\mathrm{vol}^3 H(-\dot{\varepsilon}_\mathrm{vol}) \ge 0\f$,
   * both terms non-negative -- so it can only remove energy from the system, never add it.
   *
   * The quadratic term is restricted to compression because in expansion it would resist the
   * opening of a crack or a void, which is the one place where an artificial stress is least
   * welcome. The linear term is not so restricted: damping only one sign of the volumetric rate
   * would rectify the ringing it is meant to remove into a net volumetric drift.
   *
   * @note \f$L_e\f$ must be the same characteristic length the stable time increment is computed
   * from. Both scale the term to the highest frequency the mesh can carry, and if they disagree the
   * damping is no longer tied to the mode it is meant to damp.
   */
  inline double viscousStress( double              volumetricStrainRate,
                               double              density,
                               double              waveSpeed,
                               double              characteristicElementLength,
                               const Coefficients& coefficients )
  {
    const double rateOfVolume = volumetricStrainRate;

    double stress = coefficients.linear * density * waveSpeed * characteristicElementLength * rateOfVolume;

    if ( coefficients.quadratic != 0.0 && rateOfVolume < 0.0 ) {
      const double lengthTimesB2 = coefficients.quadratic * characteristicElementLength;
      stress -= density * lengthTimesB2 * lengthTimesB2 * rateOfVolume * rateOfVolume;
    }

    return stress;
  }

  /**
   * @brief The artificial viscous stress from a strain INCREMENT rather than a rate.
   *
   * @param volumetricStrainIncrement Trace of the strain increment of the current time increment.
   * @param timeIncrement The time increment \f$\Delta t\f$ it was taken over.
   * @param density Current mass density \f$\rho\f$.
   * @param waveSpeed Dilatational wave speed \f$c_d\f$ of the material.
   * @param characteristicElementLength The element's smallest physical extent \f$L_e\f$.
   * @param coefficients The two dimensionless coefficients.
   * @return The scalar \f$\sigma_\mathrm{bv}\f$ to add to the normal stress components.
   *
   * @details Convenience overload for an explicit element, which has the strain increment at hand
   * and not the rate. A non-positive \f$\Delta t\f$ yields exactly zero rather than an infinity:
   * an explicit solver legitimately calls its elements with \f$\Delta t = 0\f$ to prime the
   * internal force at the start of a step, and a viscous stress is not defined there.
   */
  inline double viscousStressFromIncrement( double              volumetricStrainIncrement,
                                            double              timeIncrement,
                                            double              density,
                                            double              waveSpeed,
                                            double              characteristicElementLength,
                                            const Coefficients& coefficients )
  {
    if ( timeIncrement <= 0.0 )
      return 0.0;

    return viscousStress( volumetricStrainIncrement / timeIncrement,
                          density,
                          waveSpeed,
                          characteristicElementLength,
                          coefficients );
  }

} // namespace Marmot::FiniteElement::BulkViscosity
