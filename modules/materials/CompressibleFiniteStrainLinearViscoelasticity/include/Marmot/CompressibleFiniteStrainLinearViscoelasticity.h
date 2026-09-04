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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <map>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::CompressibleFiniteStrainLinearViscoelasticity
   * @brief Finite-strain compressible hyper-viscoelastic material model.
   *
   * @details The model combines a selectable hyperelastic base potential with a generalized Maxwell
   * viscoelastic evolution in the second Piola-Kirchhoff stress space.
   *
   * @par Material parameters
   * - @b baseModel - hyperelastic base identifier (@c NeoHooke, @c Yeoh, @c MooneyRivlin, @c PenceGouNeoHooke)
   * - @b onlyShearCreep - flag to restrict viscoelastic evolution to the deviatoric part
   * - @b elasticProperties - coefficients required by the selected base model
   * - @b n_Maxwell - number of Maxwell elements
   * - @b tau[i], beta[i] (i = 1..n_Maxwell) - Maxwell retardation times and relative weights
   * - @b rho - density (optional; read from the last material property entry)
   *
   * @par State variables
   * - @c S0_old - previous step stress-like internal variable used in viscoelastic update
   * - @c creepStateVars - internal variables for the generalized Maxwell chain
   */
  class CompressibleFiniteStrainLinearViscoelasticity : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /**
     * @brief Construct a CompressibleFiniteStrainLinearViscoelasticity material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    CompressibleFiniteStrainLinearViscoelasticity( const double* materialProperties,
                                                   int           nMaterialProperties,
                                                   int           materialLabel );

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
     *   - @c t - old time.
     *   - @c dT - time increment.
     *
     * Template parameter @c <3> indicates 3D.
     */
    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& ) const override;

    /**
     * @brief Get the material density.
     * @param[in] stateVars Pointer to state variables (unused).
     * @return Density from the last entry in @c materialProperties.
     */
    double getDensity( const double* stateVars ) const override;

  protected:
    /// @brief Hyperelastic base model type.
    enum HyperelasticBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, PenceGouNeoHooke = 3 };

    /// @brief Mapping from hyperelastic base model to required number of elastic properties.
    const std::map< HyperelasticBase, int > nElasticPropertiesMap = {
      { NeoHooke, 2 },
      { Yeoh, 4 },
      { MooneyRivlin, 3 },
      { PenceGouNeoHooke, 2 },
    };

    /// @brief Selected hyperelastic base model.
    const HyperelasticBase hyperelasticBase;

    /// @brief Flag indicating whether viscoelasticity acts on deviatoric stresses only.
    const double onlyShearCreep;

    /// @brief Elastic coefficients of the selected hyperelastic base model.
    const Eigen::Map< const Eigen::VectorXd > elasticProperties;

    /// @brief Generalized Maxwell model parameters.
    const ContinuumMechanics::Viscoelasticity::FiniteStrain::MaxwellProperties maxwellProperties;

    /// @brief Initial compliance tensor used for viscoelastic stress update.
    TensorUtility::FastorTensors::StandardTensors::Tensor3333d initialCompliance;

    /// @brief Define the layout of persistent state variables.
    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }

    /**
     * @brief Compute strain energy density and derivatives with respect to @f$\mathbf{C}@f$.
     * @param[in] C Right Cauchy-Green deformation tensor.
     * @return Tuple containing energy density and first, second and third derivatives.
     */
    std::tuple< double,
                TensorUtility::FastorTensors::StandardTensors::Tensor33d,
                TensorUtility::FastorTensors::StandardTensors::Tensor3333d,
                TensorUtility::FastorTensors::StandardTensors::Tensor333333d >
    computeEnergyDensityAndDerivatives( const TensorUtility::FastorTensors::StandardTensors::Tensor33d& C ) const;
  };

} // namespace Marmot::Materials
