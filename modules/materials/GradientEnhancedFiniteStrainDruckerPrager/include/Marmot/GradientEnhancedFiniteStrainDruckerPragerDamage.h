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
#include <utility>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedFiniteStrainDruckerPragerDamage
   * @brief Scalar, implicit-gradient damage of GradientEnhancedFiniteStrainDruckerPrager.
   *
   * The damage law of the finite-strain damage-plasticity models (GMCDPFiniteStrain and its gradient-enhanced
   * siblings): the LOCAL variable grows with the volumetric plastic logarithmic strain, weighted by the ductility
   * measure @f$ x_s @f$ of the compressive part of the plastic flow,
   * @f[
   *   \Delta\alpha_\mathrm{local} = \frac{\Delta\varepsilon^p_v}{x_s(R_s)},\qquad
   *   R_s = \frac{\sum_a\langle-\Delta\varepsilon^p_a\rangle}{\Delta\varepsilon^p_v},\qquad
   *   x_s = \begin{cases} 1 + A_s R_s^2 & R_s < 1\\ 1 + A_s(4\sqrt{R_s} - 3) & R_s \ge 1 \end{cases},
   * @f]
   * and the damage follows from the over-nonlocal weighting of the local and the nonlocal measure,
   * @f[
   *   \alpha_w = m\,\bar{N} + (1-m)\,\alpha_\mathrm{local},\qquad
   *   \kappa = \max_t \alpha_w,\qquad
   *   \omega = \min\bigl(1 - e^{-\kappa/\varepsilon_f},\ \omega_\mathrm{max}\bigr).
   * @f]
   * The history maximum @f$ \kappa @f$ makes the damage irreversible also when the nonlocal field decreases.
   */
  class GradientEnhancedFiniteStrainDruckerPragerDamage {
  public:
    struct ModelParameters {
      double As;                 ///< ductility parameter
      double softeningModulus;   ///< epsilon_f
      double weightingParameter; ///< m
      double maxDamage;          ///< omega_max
    };

    const ModelParameters modelParameters;

    explicit GradientEnhancedFiniteStrainDruckerPragerDamage( const ModelParameters& modelParameters )
      : modelParameters( modelParameters )
    {
    }

    /// increment of the local damage variable for principal plastic log strain increments
    double deltaAlphaLocal( const Eigen::Vector3d& deltaPlasticLogStrainPrincipal ) const;

    /// weighted measure, damage history and damage: {alphaWeighted, kappaNew, omega}
    struct Result {
      double alphaWeighted;
      double kappa;
      double omega;
    };

    Result computeDamage( double alphaLocal, double alphaNonLocal, double kappaOld ) const;

  private:
    double omega( double kappa ) const;
    double xs( double Rs ) const;
  };

} // namespace Marmot::Materials
