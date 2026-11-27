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
   * - @b baseModel - hyperelastic base identifier (@c NeoHooke, @c Yeoh, @c MooneyRivlin,
   *   @c PenceGouNeoHooke, @c ArrudaBoyce)
   * - @b onlyShearCreep - flag to restrict viscoelastic evolution to the deviatoric part
   * - @b elasticProperties - parameters of the selected compressible base potential, in the order of its
   *   interface in Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions:
   *   \f$K,G\f$ for @c NeoHooke (Coupled::StandardNeoHooke); \f$C_1,C_2,C_3,D_1\f$ for @c Yeoh
   *   (YeohPotential); \f$C_1,C_2,D_1\f$ for @c MooneyRivlin (MooneyRivlinPotential); \f$K,G\f$ for
   *   @c PenceGouNeoHooke (PenceGouPotentialB); \f$\mu,\lambda_L,D_1\f$ for @c ArrudaBoyce
   *   (ArrudaBoyce8ChainPotential). Yeoh, MooneyRivlin and ArrudaBoyce share the volumetric term
   *   \f$\frac{1}{D_1}\left(\frac{J^2-1}{2}-\ln J\right)\f$, i.e. a bulk modulus \f$2/D_1\f$.
   * - @b n_Maxwell - number of Maxwell elements
   * - @b tau[i], beta[i] (i = 1..n_Maxwell) - Maxwell retardation times and relative weights
   * - @b rho - density (optional; read from the last material property entry)
   *
   * @par State variables
   * - @c S0_old - previous step stress-like internal variable used in viscoelastic update
   * - @c creepStateVars - internal variables for the generalized Maxwell chain
   */
  class CompressibleFiniteStrainLinearViscoelasticity : public Marmot::Materials::MaterialFiniteStrain {
  public:
    using Marmot::Materials::MaterialFiniteStrain::MaterialFiniteStrain;

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
    enum HyperelasticBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, PenceGouNeoHooke = 3, ArrudaBoyce = 4 };

    /// @brief Mapping from hyperelastic base model to required number of elastic properties.
    const std::map< HyperelasticBase, int > nElasticPropertiesMap = {
      { NeoHooke, 2 },
      { Yeoh, 4 },
      { MooneyRivlin, 3 },
      { PenceGouNeoHooke, 2 },
      { ArrudaBoyce, 3 },
    };

    /// @brief Selected hyperelastic base model.
    const HyperelasticBase hyperelasticBase;

    /// @brief Flag indicating whether viscoelasticity acts on deviatoric stresses only.
    const double onlyShearCreep;

    /// @brief Elastic coefficients of the selected hyperelastic base model.
    const Eigen::Map< const Eigen::VectorXd > elasticProperties;

    /// @brief Generalized Maxwell model parameters.
    const ContinuumMechanics::Viscoelasticity::MaxwellProperties maxwellProperties;

    /// @brief Initial compliance tensor used for viscoelastic stress update.
    Marmot::Math::Tensors::FastorTensors::Tensor3333d initialCompliance;

    /// @brief Define the layout of persistent state variables.
    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }

    /**
     * @brief Initial (infinitesimal) shear modulus of the selected hyperelastic base.
     * @details @f$G_0 = 2\left(\partial\Psi/\partial\bar I_1 + \partial\Psi/\partial\bar I_2\right)@f$
     * at @f$\mathbf{C}=\mathbf{I}@f$: @f$G@f$ for @c NeoHooke and @c PenceGouNeoHooke, @f$2C_1@f$ for
     * @c Yeoh, @f$2(C_1+C_2)@f$ for @c MooneyRivlin and
     * @f$\mu\,\frac{3\lambda_L^2-1}{3(\lambda_L^2-1)}@f$ for @c ArrudaBoyce. With @c onlyShearCreep the
     * initial compliance is @f$\frac{1}{2G_0}\mathbb{I}@f$, otherwise
     * @f$\frac{1}{2G_0}\mathbb{I}^{\mathrm{dev}} + \frac{1}{9K_0}\mathbf{I}\otimes\mathbf{I}@f$.
     * @return Initial shear modulus @f$G_0@f$.
     */
    double initialShearModulus() const;

    /**
     * @brief Initial (infinitesimal) bulk modulus of the selected hyperelastic base.
     * @details @f$K_0 = \mathrm{d}^2U/\mathrm{d}J^2@f$ at @f$J=1@f$: @f$K@f$ for @c NeoHooke and
     * @c PenceGouNeoHooke, @f$2/D_1@f$ for @c Yeoh, @c MooneyRivlin and @c ArrudaBoyce.
     * @return Initial bulk modulus @f$K_0@f$.
     */
    double initialBulkModulus() const;

    /**
     * @brief Compute strain energy density and derivatives with respect to @f$\mathbf{C}@f$.
     * @param[in] C Right Cauchy-Green deformation tensor.
     * @return Tuple containing energy density and first, second and third derivatives.
     */
    std::tuple< double,
                Marmot::Math::Tensors::FastorTensors::Tensor33d,
                Marmot::Math::Tensors::FastorTensors::Tensor3333d,
                Marmot::Math::Tensors::FastorTensors::Tensor333333d >
    computeEnergyDensityAndDerivatives( const Marmot::Math::Tensors::FastorTensors::Tensor33d& C ) const;
  };

} // namespace Marmot::Materials
