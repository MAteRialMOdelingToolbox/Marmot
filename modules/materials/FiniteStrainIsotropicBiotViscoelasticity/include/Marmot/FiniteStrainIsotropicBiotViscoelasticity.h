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
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainAD.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::FiniteStrainIsotropicBiotViscoelasticity
   * @brief Finite-strain isotropic Biot-viscoelastic material model.
   *
   * @details The model combines isotropic Biot-type hyperelasticity with a generalized Maxwell evolution
   * in Biot stress space.
   *
   * @par Material parameters
   * - @b K - bulk modulus
   * - @b G - shear modulus
   * - @b n_Maxwell - number of Maxwell elements
   * - @b tau[i], beta[i] (i = 1..n_Maxwell) - Maxwell retardation times and relative weights
   * - @b rho - density (optional; read from the last material property entry)
   *
   * @par State variables
   * - @c S0_old - Biot stress from the previous increment
   * - @c creepStateVars - internal variables for the generalized Maxwell chain
   */
  class FiniteStrainIsotropicBiotViscoelasticity : public MarmotMaterialFiniteStrainAD {
  public:
    using MarmotMaterialFiniteStrainAD::MarmotMaterialFiniteStrainAD;

    /**
     * @brief Construct a FiniteStrainIsotropicBiotViscoelasticity material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    FiniteStrainIsotropicBiotViscoelasticity( const double* materialProperties,
                                              int           nMaterialProperties,
                                              int           materialLabel );

    /**
     * @brief Compute Kirchhoff stress using automatic differentiation.
     *
     * @param[in,out] response Constitutive response including stress, energy density and state variables.
     * @param[in] deformation Deformation information with deformation gradient @f$\mathbf{F}@f$.
     * @param[in] timeIncrement Time increment information.
     */
    void computeStressAD( ConstitutiveResponseAD< 3 >&,
                          const DeformationAD< 3 >&,
                          const TimeIncrement& ) const override;

    /**
     * @brief Get the material density.
     * @param[in] stateVars Pointer to state variables (unused).
     * @return Density from the last entry in @c materialProperties.
     */
    double getDensity( const double* stateVars ) const override;

  protected:
    /// @brief Bulk modulus.
    const double K;
    /// @brief Shear modulus.
    const double G;

    /// @brief Generalized Maxwell model parameters.
    const ContinuumMechanics::Viscoelasticity::FiniteStrain::MaxwellProperties maxwellProperties;

    /// @brief Initial compliance tensor used for viscoelastic stress update.
    const TensorUtility::FastorTensors::StandardTensors::Tensor3333t< autodiff::dual > initialCompliance;

    /// @brief Define the layout of persistent state variables.
    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }
  };

} // namespace Marmot::Materials
