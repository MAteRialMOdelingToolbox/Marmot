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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <string>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::CompressibleNeoHooke
   * @brief Compressible Neo-Hookean hyperelastic material model (Pence–Gou potential, variant B).
   *
   * @par Material parameters
   * - @b K - bulk modulus
   * - @b G - shear modulus
   *
   * @par State variables
   * - No state variables required.
   *
   * @ingroup materials_hyperelastic
   */

  class CompressibleNeoHooke : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /**
     * @brief Construct a CompressibleNeoHooke material.
     * @param materialProperties Expects @c K at index 0 and @c G at index 1.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */

    CompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    static constexpr int nStateVarsRequired = 0; /**< Number of required state variables (none here). */

    /**
     * @brief Compute the Kirchhoff stress and the algorithmic tangent for the current step.
     *
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$.
     *   - @c elasticEnergyDensity - elastic energy density  @f$\psi@f$.
     *   - @c rho - density (unused here).
     * @param[in,out] tangents
     *   - @c dTau_dF - algorithmic tangent @f$\partial\boldsymbol{\tau}/\partial\boldsymbol{F}@f$.
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement
     *   - @c t - old (pseudo-)time.
     *   - @c dT - (pseudo-)time increment.
     *
     * Template parameter @c <3> indicates 3D.
     */

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& );
    /** @brief Number of required state variables.
     *  @return Always 0 for this model.
     */
    int getNumberOfRequiredStateVars() { return this->nStateVarsRequired; }

    /** @brief Bind external state storage (unused for this model; required for the interface).
     *  @param stateVars Pointer to a contiguous array provided by the caller for internal state.
     *  @param nStateVars Number of entries in that array.
     */
    void assignStateVars( double* stateVars, int nStateVars )
    {
      this->stateVars  = stateVars;
      this->nStateVars = nStateVars;
    };

    /**
     * @brief Access a named state quantity (no states here).
     * @param result Name of the state to view.
     * @return Always an empty StateView since no states are used.
     */

    StateView getStateView( const std::string& result );
  };

} // namespace Marmot::Materials
