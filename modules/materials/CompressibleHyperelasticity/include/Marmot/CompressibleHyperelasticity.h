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
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <variant>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::CompressibleHyperelasticity
   * @brief Compressible hyperelastic material model with a decoupled potential
   * \f$\Psi(\boldsymbol C) = \Psi_{\rm iso}(\bar I_1, \bar I_2) + U(J)\f$, whose isochoric and volumetric parts are
   * selected by the first two material properties.
   *
   * @par Material properties
   * - @b isochoricBase - isochoric part, see #IsochoricBase
   * - @b volumetricBase - volumetric part, see #VolumetricBase
   * - parameters of the isochoric part, in the order of its type in
   *   Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric:
   *   \f$\mu\f$ for @c NeoHooke; \f$C_{10}, C_{20}, C_{30}\f$ for @c Yeoh; \f$C_{10}, C_{01}\f$ for @c MooneyRivlin;
   *   \f$\mu, \lambda_L\f$ for @c ArrudaBoyce8Chain
   * - parameters of the volumetric part, in the order of its type in
   *   Marmot::ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric:
   *   \f$\kappa, \theta\f$ for @c Ogden; \f$\kappa, \alpha, \beta\f$ for @c DollSchweizerhof; \f$\kappa\f$ for all
   *   others
   * - @b rho - density (optional, only read by getDensity())
   *
   * The compressible Neo-Hooke model of Pence & Gou (variant B) with bulk modulus \f$K\f$ and shear modulus
   * \f$G\f$ is the pair @c NeoHooke and @c DollSchweizerhof with the properties
   * \f$[0, 6, G, K, 1, 3]\f$.
   *
   * @par State variables
   * - No state variables required.
   */
  class CompressibleHyperelasticity : public Marmot::Materials::MaterialFiniteStrain {
  public:
    using Marmot::Materials::MaterialFiniteStrain::MaterialFiniteStrain;

    /// @brief Isochoric part of the potential, selected by the first material property.
    enum IsochoricBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, ArrudaBoyce8Chain = 3 };

    /// @brief Volumetric part of the potential, selected by the second material property; numbered as
    /// \f$U_1\f$ to \f$U_8\f$ in Doll & Schweizerhof (2000), https://doi.org/10.1115/1.321146.
    enum VolumetricBase {
      SussmanBathe                = 1,
      Penalty                     = 2,
      SimoTaylor                  = 3,
      Ogden                       = 4,
      LiuHofstetterMang           = 5,
      DollSchweizerhof            = 6,
      DollSchweizerhofExponential = 7,
      DollSchweizerhofLogarithmic = 8
    };

    /**
     * @brief Construct a CompressibleHyperelasticity material.
     * @param materialProperties Isochoric base, volumetric base, the parameters of the isochoric part, the
     * parameters of the volumetric part and, optionally, the density.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     * @throws std::invalid_argument for an unknown base, too few material properties or an Arruda-Boyce locking
     * stretch \f$\lambda_L \le 1\f$.
     */
    CompressibleHyperelasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

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
                        const TimeIncrement& ) const override;

    /**
     * @brief Get material density.
     * @param stateVars State variables (unused).
     * @return Density value, the material property following the parameters of the volumetric part.
     */
    double getDensity( const double* stateVars ) const override;

  protected:
    /// @brief Scalar function of the isochoric part.
    using IsochoricPart = std::variant<
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric::NeoHooke,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric::Yeoh,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric::MooneyRivlin,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Isochoric::ArrudaBoyce8Chain >;

    /// @brief Scalar function of the volumetric part.
    using VolumetricPart = std::variant<
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::SussmanBathe,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::Penalty,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::SimoTaylor,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::Ogden,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::LiuHofstetterMang,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::DollSchweizerhof,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::DollSchweizerhofExponential,
      ContinuumMechanics::Elasticity::EnergyDensityFunctions::Volumetric::DollSchweizerhofLogarithmic >;

    /// @brief Isochoric part with its parameters.
    IsochoricPart isochoricPart;

    /// @brief Volumetric part with its parameters.
    VolumetricPart volumetricPart;

    /// @brief Index of the density in the material properties.
    int indexDensity = 0;
  };

} // namespace Marmot::Materials
