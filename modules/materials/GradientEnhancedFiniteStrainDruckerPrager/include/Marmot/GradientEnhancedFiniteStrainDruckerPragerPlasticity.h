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
 * Thomas Mader thomas.mader@boku.ac.at
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
#include <Eigen/Core>
#include <exception>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedFiniteStrainDruckerPragerPlasticity
   * @brief Drucker-Prager return mapping of GradientEnhancedFiniteStrainDruckerPrager.
   *
   * Works in the PRINCIPAL elastic logarithmic strains @f$ \varepsilon^e_a = \ln\lambda^e_a @f$ of the
   * elastic left stretch. For the isotropic compressible neo-Hookean potential of CompressibleNeoHooke,
   * @f[
   *   \psi = \frac{K}{8}\bigl(J^2 + J^{-2} - 2\bigr) + \frac{G}{2}\bigl(I_1 J^{-2/3} - 3\bigr),
   * @f]
   * the Mandel stress is coaxial with the elastic stretch, and its principal values are
   * @f[
   *   \Sigma_a = \frac{\partial\psi}{\partial\varepsilon^e_a}
   *            = \frac{K}{2}\sinh(2\theta) + G\Bigl(e^{2 e_a} - \tfrac13\sum_b e^{2 e_b}\Bigr),
   *   \qquad \theta = \textstyle\sum_a\varepsilon^e_a,\quad e_a = \varepsilon^e_a - \theta/3 ,
   * @f]
   * which coincide with the principal Kirchhoff stresses. With the exponential map of the plastic flow,
   * the multiplicative return map is then exactly the small-strain one in these variables:
   * @f[
   *   \varepsilon^e_a = \varepsilon^{e,\mathrm{trial}}_a - \Delta\lambda\,\frac{\partial g}{\partial\Sigma_a},\qquad
   *   f(\boldsymbol{\Sigma}, \alpha) = \sqrt{J_2} + \eta\,p - \xi\,c(\alpha) = 0,\qquad
   *   \alpha = \alpha_n + \xi\,\Delta\lambda ,
   * @f]
   * with @f$ p = \tfrac13\mathrm{tr}\,\boldsymbol\Sigma @f$ (tension positive), the plastic potential
   * @f$ g = \sqrt{J_2} + \bar\eta\,p @f$ and linear hardening of the cohesion, @f$ c = c_0 + H\alpha @f$.
   * When the deviatoric stress would reverse, the state is returned to the APEX of the cone instead
   * (de Souza Neto, Peric & Owen, Computational Methods for Plasticity, Sec. 8.3), where
   * @f$ \Delta\varepsilon^p_v @f$ is the unknown and @f$ \alpha = \alpha_n + (\xi/\bar\eta)\,\Delta\varepsilon^p_v @f$.
   */
  class GradientEnhancedFiniteStrainDruckerPragerPlasticity {
  public:
    struct ReturnMappingFailedException : std::exception {};

    struct MaterialParameters {
      double K; ///< bulk modulus
      double G; ///< shear modulus
    };

    struct ModelParameters {
      double eta;    ///< friction parameter of the yield function
      double xi;     ///< cohesion parameter of the yield function
      double etaBar; ///< dilatancy parameter of the plastic potential
      double c0;     ///< initial cohesion
      double H;      ///< linear hardening modulus of the cohesion
    };

    struct StressState {
      Eigen::Vector3d mandelPrincipal; ///< principal Mandel (= Kirchhoff) stresses
    };

    struct MaterialState {
      Eigen::Vector3d elasticLogStrainPrincipal;
      double          alphaP; ///< hardening variable
    };

    enum class ReturnMode { Elastic, Cone, Apex };

    struct ReturnMapResult {
      StressState     newStressState;
      MaterialState   newMaterialState;
      Eigen::Vector3d deltaPlasticLogStrainPrincipal; ///< trial minus new elastic log strains
      ReturnMode      mode;
    };

    GradientEnhancedFiniteStrainDruckerPragerPlasticity( const MaterialParameters&, const ModelParameters& );

    StressState computeStressState( const MaterialState& ) const;

    /// @f$ \partial\Sigma_a / \partial\varepsilon^e_b @f$
    Eigen::Matrix3d dStress_dElasticLogStrain( const MaterialState& ) const;

    double yieldFunction( const StressState&, double alphaP ) const;

    bool checkIfYielding( const StressState&, const MaterialState& ) const;

    ReturnMapResult performReturnMapping( const MaterialState& trialState ) const;

  private:
    const MaterialParameters materialParameters;
    const ModelParameters    modelParameters;

    double cohesion( double alphaP ) const { return modelParameters.c0 + modelParameters.H * alphaP; }

    bool returnToCone( const MaterialState& trialState, ReturnMapResult& result ) const;
    bool returnToApex( const MaterialState& trialState, ReturnMapResult& result ) const;
  };

} // namespace Marmot::Materials
