/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include <array>

namespace Marmot::Materials {

  /**
   * @brief Gradient-enhanced (implicit gradient) version of RankineDamagePlasticity: multi-surface Rankine plasticity
   * in the effective stress with isotropic damage driven by a nonlocal variable, following the damage-plastic model of
   * Grassl and Jirásek (2006, IJSS 43, 7166-7196) and the structure of GCDP.
   *
   * @details The effective stress \f$ \bar{\boldsymbol{\sigma}} \f$ is bounded by the three Rankine surfaces of the
   * perfectly plastic RankineSofteningPlasticity, and the local damage driving variable is the volumetric plastic
   * strain \f$ \kappa_\mathrm{d} \f$ with \f$ \Delta\kappa_\mathrm{d} = \sum_A \Delta\lambda_A \f$, as in
   * RankineDamagePlasticity. It is the source of the Helmholtz equation solved by the element,
   * \f[
   *   \bar\kappa_\mathrm{d} - l^2\,\Delta\bar\kappa_\mathrm{d} = \kappa_\mathrm{d}, \qquad
   *   \nabla\bar\kappa_\mathrm{d}\cdot\boldsymbol{n} = 0 \text{ on the boundary},
   * \f]
   * and the damage follows from the over-nonlocal combination (\f$ m = 1 \f$: standard nonlocal, \f$ m > 1 \f$:
   * over-nonlocal, \f$ m = 0 \f$: local)
   * \f[
   *   \hat\kappa = m\,\bar\kappa_\mathrm{d} + (1-m)\,\kappa_\mathrm{d}, \qquad
   *   \omega = 1 - \exp(-\hat\kappa_\mathrm{max}/\varepsilon_\mathrm{f}), \qquad
   *   \hat\kappa_\mathrm{max} = \max_{\tau \le t} \hat\kappa(\tau), \qquad
   *   \boldsymbol{\sigma} = (1-\omega)\,\bar{\boldsymbol{\sigma}}.
   * \f]
   * The history maximum keeps the damage irreversible also for \f$ m > 1 \f$, where \f$ \hat\kappa \f$ may decrease.
   * Since the effective stress does not depend on the damage, \f$ \kappa_\mathrm{d} \f$ does not depend on
   * \f$ \bar\kappa_\mathrm{d} \f$. Free energy \f$ (1-\omega)\,\bar\psi_\mathrm{el} \f$ and dissipation increment
   * \f$ (1-\omega)\,\bar{\boldsymbol{\sigma}}:\Delta\boldsymbol{\varepsilon}^\mathrm{p} + \bar\psi_\mathrm{el}\,
   * \Delta\omega \f$ as in RankineDamagePlasticity.
   *
   * Returned to the element: \f$ \kappa_\mathrm{d} \f$ (KLocal), \f$ c = l^2 \f$, and the tangents
   * \f$ \partial\boldsymbol{\sigma}/\partial\Delta\boldsymbol{\varepsilon} = (1-\omega)\,\bar{\mathbb{C}}^\mathrm{ep} -
   * (1-m)\,\omega'\,\bar{\boldsymbol{\sigma}} \otimes \partial\kappa_\mathrm{d}/\partial\Delta\boldsymbol{\varepsilon}
   * \f$, \f$ \partial\boldsymbol{\sigma}/\partial\bar\kappa_\mathrm{d} = -m\,\omega'\,\bar{\boldsymbol{\sigma}} \f$,
   * \f$ \partial\kappa_\mathrm{d}/\partial\Delta\boldsymbol{\varepsilon} = \boldsymbol{1} -
   * \bar{\mathbb{C}}^{\mathrm{ep},T}:\mathbb{C}^{-1}:\boldsymbol{1} \f$ and \f$ \partial\kappa_\mathrm{d}/\partial
   * \bar\kappa_\mathrm{d} = 0 \f$, with \f$ \omega' = \exp(-\hat\kappa/\varepsilon_\mathrm{f})/\varepsilon_\mathrm{f}
   * \f$ while \f$ \hat\kappa \f$ exceeds its history maximum, and zero otherwise.
   *
   * Material properties: \f$ E, \nu, f_\mathrm{t}, \varepsilon_\mathrm{f}, l, m \f$, optionally the density and the
   * nonlocal viscosity (default zero, for explicit dynamics only).
   */
  class RankineGradientDamagePlasticity : public MaterialGradientEnhancedHypoElastic< 1 > {

  public:
    using Base = MaterialGradientEnhancedHypoElastic< 1 >;
    using typename Base::increment;
    using typename Base::response;
    using typename Base::tangents;

    RankineGradientDamagePlasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

    RankineGradientDamagePlasticity( const RankineGradientDamagePlasticity& )            = delete;
    RankineGradientDamagePlasticity& operator=( const RankineGradientDamagePlasticity& ) = delete;

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    /// the stress update of computeStress, without the check of its results
    void updateStress( response& res, tangents& tan, const increment& inc ) const;

    /// properties of the perfectly plastic RankineSofteningPlasticity for the effective stress
    const std::array< double, 6 >    effectiveProperties;
    const RankineSofteningPlasticity effectiveMaterial;
  };

} // namespace Marmot::Materials
