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
   * @class Marmot::Materials::FiniteStrainOrthotropicBiotViscoelasticity
   * @brief Finite-strain orthotropic Biot-viscoelastic material model.
   *
   * @details The model combines orthotropic Biot-type hyperelasticity with a generalized Maxwell evolution
   * in Biot stress space.
   *
   * @par Material parameters
   * - @b E1, @b E2, @b E3 - Young's moduli in principal material directions
   * - @b nu12, @b nu13, @b nu23 - Poisson ratios
   * - @b G12, @b G13, @b G23 - shear moduli
   * - @b n_Maxwell - number of Maxwell elements
   * - @b tau[i], beta[i] (i = 1..n_Maxwell) - Maxwell retardation times and relative weights
   * - @b rho - density (optional; read from the last material property entry)
   *
   * @par State variables
   * - @c S0_old - Biot stress from the previous increment
   * - @c creepStateVars - internal variables for the generalized Maxwell chain
   */
  class FiniteStrainOrthotropicBiotViscoelasticity : public MarmotMaterialFiniteStrainAD {
  public:
    using MarmotMaterialFiniteStrainAD::MarmotMaterialFiniteStrainAD;

    /**
     * @brief Construct a FiniteStrainOrthotropicBiotViscoelasticity material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    FiniteStrainOrthotropicBiotViscoelasticity( const double* materialProperties,
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
    /// @brief Young's modulus in principal material direction 1.
    const double E1;
    /// @brief Young's modulus in principal material direction 2.
    const double E2;
    /// @brief Young's modulus in principal material direction 3.
    const double E3;

    /// @brief Poisson ratio relating directions 1 and 2.
    const double nu12;
    /// @brief Poisson ratio relating directions 1 and 3.
    const double nu13;
    /// @brief Poisson ratio relating directions 2 and 3.
    const double nu23;

    /// @brief Shear modulus in the 1-2 plane.
    const double G12;
    /// @brief Shear modulus in the 1-3 plane.
    const double G13;
    /// @brief Shear modulus in the 2-3 plane.
    const double G23;

    /// @brief Generalized Maxwell model parameters.
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellProperties;

    /// @brief Constant derivative of Biot stress with respect to right stretch tensor.
    const FastorStandardTensors::Tensor3333t< autodiff::dual > dBiotStress_dU;

    /// @brief Define the layout of persistent state variables.
    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }
  };

} // namespace Marmot::Materials
