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
   * \class Marmot::Materials::CompressibleNeoHooke
   * \brief Compressible Neo-Hookean hyperelastic material using the Pence–Gou potential (variant B).
   *
   * Computes Kirchhoff stress \f$\boldsymbol{\tau}\f$ and the algorithmic tangent
   * \f$\partial\boldsymbol{\tau}/\partial\mathbf{F}\f$ from the deformation gradient \f$\mathbf{F}\f$. The strain
   * energy is evaluated via `EnergyDensityFunctions::PenceGouPotentialB(K,G)`.
   *
   * \par Material parameters (indices in #materialProperties)
   * - \c K (#materialProperties[0]) — bulk modulus [Pa]
   * - \c G (#materialProperties[1]) — shear modulus [Pa]
   *
   * \par State variables
   * None. \c nStateVarsRequired = 0.
   *
   * \par Outputs set on \c response
   * - \c response.tau : Kirchhoff stress
   * - \c response.elasticEnergyDensity : strain energy density
   * - \c response.rho : set to 1.0 (density scaling not handled here)
   *
   * \ingroup materials_hyperelastic
   */

  class CompressibleNeoHooke : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /**
     * \brief Construct a CompressibleNeoHooke material.
     * \param materialProperties Array with at least 2 entries: \c K at [0], \c G at [1].
     * \param nMaterialProperties Length of \p materialProperties.
     * \param materialLabel User-defined material label (passed to base).
     */

    CompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    static constexpr int nStateVarsRequired = 0; /**< Number of required state variables (none). */

    /**
     * \brief Compute Kirchhoff stress and elastic tangent for the current configuration.
     *
     * Uses \f$\mathbf{C} = \mathbf{F}^\mathrm{T}\mathbf{F}\f$ and the Pence–Gou energy to obtain
     * \f$\boldsymbol{\tau}\f$ and \f$\partial\boldsymbol{\tau}/\partial\mathbf{F}\f$.
     *
     * \param[out] response  Filled with \c tau, \c elasticEnergyDensity, \c rho.
     * \param[out] tangents  Filled with \c dTau_dF.
     * \param[in]  deformation  Contains deformation gradient tensor F.
     * \param[in]  timeIncrement  Current time step information \c t and \c dT.
     * \note No state variables are updated (purely elastic model).
     */

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& );

    int getNumberOfRequiredStateVars() { return this->nStateVarsRequired; } /**< Always returns 0. */

    /** \brief Attach external state storage (unused for this model; required for the interface).
     *  \param stateVars Unused.
     *  \param nStateVars Unused.
     */
    void assignStateVars( double* stateVars, int nStateVars )
    {
      this->stateVars  = stateVars;
      this->nStateVars = nStateVars;
    };

    /**
     * \brief Access a named state quantity.
     * \param result Name of the state to view.
     * \return A view into the requested state.
     * \throws std::out_of_range If the state name is unknown (this model defines no states).
     */

    StateView getStateView( const std::string& result );
  };

} // namespace Marmot::Materials
