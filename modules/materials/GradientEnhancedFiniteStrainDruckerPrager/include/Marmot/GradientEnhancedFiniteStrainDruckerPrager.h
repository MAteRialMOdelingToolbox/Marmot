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
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include <utility>

namespace Marmot::Materials {

  /**
   * @class GradientEnhancedFiniteStrainDruckerPrager
   * @brief Finite-strain Drucker-Prager damage-plasticity with implicit-gradient (nonlocal) damage.
   *
   * **Kinematics.** Multiplicative, @f$ \boldsymbol{F} = \boldsymbol{F}^e\boldsymbol{F}^p @f$, with the plastic
   * deformation gradient as state. The trial elastic state @f$ \boldsymbol{F}^{e,\mathrm{trial}} =
   * \boldsymbol{F}\,{\boldsymbol{F}^p_n}^{-1} @f$ is decomposed spectrally through its left Cauchy-Green tensor,
   * the plastic flow is integrated with the exponential map, and the elastic rotation is frozen over the step
   * (as in GMCDPFiniteStrain). For isotropic elasticity and an isotropic yield function this makes the return map
   * EXACTLY the small-strain one in the principal elastic logarithmic strains, see
   * GradientEnhancedFiniteStrainDruckerPragerPlasticity; the plastic deformation gradient is updated as
   * @f$ \boldsymbol{F}^p_{n+1} =
   * {\boldsymbol{F}^{e,\mathrm{trial}}}^{-1}\exp(\Delta\boldsymbol{\varepsilon}^p)\,\boldsymbol{F}
   * @f$, so that @f$ \boldsymbol{F}\,{\boldsymbol{F}^p_{n+1}}^{-1} = \boldsymbol{V}^e_{n+1}\boldsymbol{R}^e @f$.
   *
   * **Elasticity** is that of CompressibleNeoHooke (PenceGouPotentialB) in the elastic stretch.
   *
   * **Plasticity.** Drucker-Prager yield function on the Mandel (= Kirchhoff, by isotropy) stress,
   * @f$ f = \sqrt{J_2} + \eta\,p - \xi\,(c_0 + H\alpha) @f$, non-associated flow with the dilatancy
   * @f$ \bar\eta @f$, and an apex return. The parameters are those of the Drucker-Prager cone through the OUTER
   * edges of the Mohr-Coulomb pyramid (compressive meridian),
   * @f[
   *   \eta = \frac{6\sin\phi}{\sqrt3\,(3-\sin\phi)},\qquad \xi = \frac{6\cos\phi}{\sqrt3\,(3-\sin\phi)},\qquad
   *   \bar\eta = \frac{6\sin\psi}{\sqrt3\,(3-\sin\psi)} .
   * @f]
   *
   * **Damage** is the implicit-gradient damage of the finite-strain damage-plasticity models, see
   * GradientEnhancedFiniteStrainDruckerPragerDamage: the local variable (volumetric plastic log strain
   * over the ductility measure) is the source @f$ L @f$ of the nonlocal balance, the damage follows from
   * @f$ m\bar{N} + (1-m)\alpha_\mathrm{local} @f$, and @f$ \boldsymbol{\tau} =
   * (1-\omega)\,\boldsymbol{\tau}_\mathrm{eff} @f$.
   *
   * The algorithmic tangents are computed by forward finite differences of the full state update, with the state
   * frozen at the beginning of the increment for every probe.
   *
   * Material properties:
   * | idx | symbol                  | meaning                                            |
   * |-----|-------------------------|----------------------------------------------------|
   * | 0   | @f$ K @f$               | bulk modulus                                       |
   * | 1   | @f$ G @f$               | shear modulus                                      |
   * | 2   | @f$ c_0 @f$             | cohesion                                           |
   * | 3   | @f$ \phi @f$            | friction angle [deg]                               |
   * | 4   | @f$ \psi @f$            | dilatancy angle [deg]                              |
   * | 5   | @f$ H @f$               | linear hardening modulus of the cohesion           |
   * | 6   | @f$ A_s @f$             | ductility parameter of the damage                  |
   * | 7   | @f$ \varepsilon_f @f$   | softening modulus of the damage                    |
   * | 8   | @f$ \omega_\max @f$     | maximum damage                                     |
   * | 9   | @f$ l @f$               | nonlocal radius, @f$ c = l^2 @f$                   |
   * | 10  | @f$ m @f$               | weighting of the nonlocal measure in the damage    |
   * | 11  | @f$ \rho @f$            | density in the reference configuration (optional)  |
   *
   * State variables: @c Fp (9), @c alphaP (hardening variable), @c alphaD (local damage variable),
   * @c kappa (damage history), @c omega (damage).
   */
  class GradientEnhancedFiniteStrainDruckerPrager : public MarmotMaterialGradientEnhancedFiniteStrain {

  public:
    GradientEnhancedFiniteStrainDruckerPrager( const double* materialProperties,
                                               int           nMaterialProperties,
                                               int           materialNumber );

    using MarmotMaterialGradientEnhancedFiniteStrain::computeStress;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override;

    /// Fp must be initialized to the identity (the base default zero-fill would be singular).
    void initializeYourself( double* stateVars, int nStateVars ) override;

    /// Drucker-Prager parameters {eta, xi} for a Mohr-Coulomb angle in degrees (outer cone).
    static std::pair< double, double > outerConeParameters( double angleInDegrees );

  protected:
    const double& K;
    const double& G;
    const double& c0;
    const double& frictionAngle;
    const double& dilatancyAngle;
    const double& H;
    const double& As;
    const double& softeningModulus;
    const double& maxDamage;
    const double& nonLocalRadius;
    const double& weightingParameter;

    /// the full state update for (F, Nbar), written into the state vector sv: {tau, L, elastic energy density}
    std::tuple< Fastor::Tensor< double, 3, 3 >, double, double > stressUpdate( const Fastor::Tensor< double, 3, 3 >& F,
                                                                               double  nonLocalField,
                                                                               double* sv ) const;
  };

} // namespace Marmot::Materials
