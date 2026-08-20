/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include "Marmot/MarmotMaterialFiniteStrainAD.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::Ogden
   * @brief Classical N-term Ogden hyperelastic material model, purely isochoric (incompressible).
   *
   * @par Material parameters
   * - @b nOgdenTerms - number of Ogden terms \f$N\f$ (materialProperties[0])
   * - @b mu_1 ... mu_N - Ogden moduli \f$\mu_p\f$ (materialProperties[1] ... materialProperties[N])
   * - @b alpha_1 ... alpha_N - Ogden exponents \f$\alpha_p\f$ (materialProperties[N+1] ...
   *   materialProperties[2N])
   * - @b density (optional) - materialProperties[2N+1]
   *
   * @par State variables
   * - No state variables required.
   *
   * @details Implements only the isochoric part of the classical (multi-term) Ogden potential,
   * \f[
   *   \Psi_\mathrm{iso}(\bar\lambda) = \sum_{p=1}^N \frac{\mu_p}{\alpha_p} \left(
   *   \bar\lambda_1^{\alpha_p} + \bar\lambda_2^{\alpha_p} + \bar\lambda_3^{\alpha_p} - 3 \right),
   * \f]
   * with isochoric principal stretches \f$\bar\lambda_i = J^{-1/3}\lambda_i\f$. There is
   * intentionally no volumetric energy term \f$U(J)\f$: this material resists no volumetric
   * deformation on its own and must therefore not be used with a plain displacement element. It
   * is intended exclusively for use with a mixed displacement-pressure element that enforces
   * incompressibility as an independent constraint, e.g.
   * Marmot::Elements::DisplacementPressureFiniteStrainElement.
   *
   * @ingroup materials_hyperelastic
   */
  class Ogden : public MarmotMaterialFiniteStrainAD {

  public:
    /**
     * @brief Construct an Ogden material.
     * @param materialProperties Expects @c nOgdenTerms at index 0, followed by the @c N moduli
     * @c mu_p, followed by the @c N exponents @c alpha_p, optionally followed by a density value.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    Ogden( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * @brief Compute the isochoric Kirchhoff stress with dual numbers.
     *
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$ (trace-free/deviatoric).
     *   - @c elasticEnergyDensity - isochoric elastic energy density @f$\Psi_\mathrm{iso}@f$.
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement Unused (hyperelastic, no history).
     *
     * Template parameter @c <3> indicates 3D.
     */
    void computeStressAD( ConstitutiveResponseAD< 3 >& response,
                          const DeformationAD< 3 >&    deformation,
                          const TimeIncrement&         timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    const int nOgdenTerms;
  };

} // namespace Marmot::Materials
