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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::IncompressibleMooneyRivlin
   * @brief Classical incompressible (purely isochoric) Mooney-Rivlin hyperelastic material model.
   *
   * @par Material parameters
   * - @b C1 - first Mooney-Rivlin modulus (materialProperties[0])
   * - @b C2 - second Mooney-Rivlin modulus (materialProperties[1])
   * - @b density (optional) - materialProperties[2]
   *
   * @par State variables
   * - No state variables required.
   *
   * @details Implements the standard isochoric Mooney-Rivlin potential
   * \f[
   *   \Psi_\mathrm{iso}(\bar\lambda) = C_1 \left( \bar I_1 - 3 \right) + C_2 \left( \bar I_2 - 3
   *   \right),
   * \f]
   * with isochoric principal stretches \f$\bar\lambda_i = J^{-1/3}\lambda_i\f$,
   * \f$\bar I_1 = \bar\lambda_1^2+\bar\lambda_2^2+\bar\lambda_3^2 = I_1 J^{-2/3}\f$ and
   * \f$\bar I_2 = I_2 J^{-4/3}\f$, with \f$I_1=\operatorname{tr}\boldsymbol C\f$ and
   * \f$I_2=\tfrac12\left(I_1^2-\operatorname{tr}\boldsymbol C^2\right)\f$. Although mathematically
   * equivalent to the two-term Ogden model with \f$(\mu_1,\alpha_1)=(2C_1,2)\f$ and
   * \f$(\mu_2,\alpha_2)=(-2C_2,-2)\f$ (see Marmot::Materials::Ogden), the potential is quadratic in
   * the stretches and is therefore evaluated directly from the invariants of \f$\boldsymbol C\f$
   * (EnergyDensityFunctions::Incompressible::MooneyRivlinPotential()) - no spectral decomposition
   * of \f$\boldsymbol C\f$ is performed, matching the efficiency of the compressible Mooney-Rivlin
   * potential.
   *
   * The consistent algorithmic tangent is computed analytically via
   * EnergyDensityFunctions::Incompressible::SecondOrderDerived::MooneyRivlinPotential() - no
   * automatic differentiation is used, matching Marmot::Materials::CompressibleNeoHooke's approach.
   *
   * There is no volumetric energy term \f$U(J)\f$: this material resists no volumetric deformation
   * on its own and must therefore not be used with a plain displacement element. It is intended
   * exclusively for use with a mixed displacement-pressure element that enforces incompressibility
   * as an independent constraint, e.g. Marmot::Elements::DisplacementPressureFiniteStrainElement.
   *
   * @ingroup materials_hyperelastic
   */
  class IncompressibleMooneyRivlin : public MarmotMaterialFiniteStrain {

  public:
    /**
     * @brief Construct an IncompressibleMooneyRivlin material.
     * @param materialProperties Expects @c C1 at index 0 and @c C2 at index 1.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    IncompressibleMooneyRivlin( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * @brief Compute the isochoric Kirchhoff stress and the algorithmic tangent for the current step.
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$ (trace-free/deviatoric).
     *   - @c elasticEnergyDensity - isochoric elastic energy density @f$\Psi_\mathrm{iso}@f$.
     * @param[in,out] tangents
     *   - @c dTau_dF - algorithmic tangent @f$\partial\boldsymbol{\tau}/\partial\boldsymbol{F}@f$.
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement Unused (hyperelastic, no history).
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override
    {
      if ( this->nMaterialProperties < 3 )
        throw std::runtime_error(
          std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 3." ) );
      return this->materialProperties[2];
    }
  };

} // namespace Marmot::Materials
