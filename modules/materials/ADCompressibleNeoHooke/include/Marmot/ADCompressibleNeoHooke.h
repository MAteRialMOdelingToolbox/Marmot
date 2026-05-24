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
#include "Marmot/MarmotMaterialFiniteStrainAD.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::ADCompressibleNeoHooke
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

  class ADCompressibleNeoHooke : public MarmotMaterialFiniteStrainAD {
  public:
    using MarmotMaterialFiniteStrainAD::MarmotMaterialFiniteStrainAD;

    /**
     * @brief Construct a ADCompressibleNeoHooke material.
     * @param materialProperties Expects @c K at index 0 and @c G at index 1.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    ADCompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * @brief Compute the Kirchhoff stress with dual numbers.
     *
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$.
     *   - @c elasticEnergyDensity - elastic energy density  @f$\psi@f$.
     *   - @c rho - density (unused here).
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement
     *   - @c t - old (pseudo-)time.
     *   - @c dT - (pseudo-)time increment.
     *
     * Template parameter @c <3> indicates 3D.
     */
    void computeStressAD( ConstitutiveResponseAD< 3 >&,
                          const DeformationAD< 3 >&,
                          const TimeIncrement& ) const override;

    double getDensity( const double* stateVars ) const override
    {
      if ( this->nMaterialProperties < 3 ) {
        throw std::runtime_error(
          std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 3." ) );
      }
      return this->materialProperties[2];
    }
  };

} // namespace Marmot::Materials
