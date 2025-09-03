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

/** \file CompressibleNeoHooke.h
 *  \brief Compressible Neo-Hookean hyperelastic material (finite strain).
 *  \details Purely elastic model: computes Kirchhoff stress \f$\tau\f$ and consistent
 *  algorithmic moduli from a compressible energy density (Pence–Gou potential B).
 *  Material parameters are read from \c materialProperties as:
 *   - [0] \c K — bulk modulus [MPa]
 *   - [1] \c G — shear modulus [MPa]
 *  No history variables are used by this model.
 */

#pragma once
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <string>

namespace Marmot::Materials {

  /** \class Marmot::Materials::CompressibleNeoHooke
   *  \brief Compressible Neo-Hookean hyperelastic material.
   *  \details Implements stress and tangent via derivatives of the energy w.r.t.
   *  the right Cauchy–Green tensor \f$C=F^\mathrm{T}F\f$. State vector is unused.
   */

  class CompressibleNeoHooke : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /** \brief Construct the compressible Neo-Hooke material.
     *  \param materialProperties Expected order: \c { K, G } (both in MPa).
     *  \param nMaterialProperties Number of entries in \c materialProperties (>=2).
     *  \param materialLabel Identifier forwarded to the base class.
     */

    CompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    static constexpr int nStateVarsRequired = 0;

    /** \brief Compute Kirchhoff stress and consistent algorithmic moduli.
     *  \param[out] response  Fills \c tau (Kirchhoff stress), \c rho (if used by framework),
     *                        and \c elasticEnergyDensity.
     *  \param[out] tangents  Consistent \f$\partial \tau / \partial F\f$.
     *  \param[in]  deformation Uses \c deformation.F.
     *  \param[in]  timeIncrement Provided for interface compatibility (not used here).
     */

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& );

    /** \brief Number of required state variables (always 0 for this elastic model). */

    int getNumberOfRequiredStateVars() { return this->nStateVarsRequired; }

    /** \brief Assign external state storage (kept for interface compatibility; not used). */

    void assignStateVars( double* stateVars, int nStateVars )
    {
      this->stateVars  = stateVars;
      this->nStateVars = nStateVars;
    };

    /** \brief Return a named state view (none are defined for this model).
     *  \throws std::out_of_range if \p stateName is unknown.
     */

    StateView getStateView( const std::string& result );
  };

} // namespace Marmot::Materials
