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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/RankineSofteningPlasticity.h"
#include <array>

namespace Marmot::Materials {

  /**
   * @brief Multi-surface Rankine plasticity in the effective stress combined with isotropic damage, following the
   * damage-plastic model of Grassl and Jirásek (2006, IJSS 43, 7166-7196).
   *
   * @details The effective stress \f$ \bar{\boldsymbol{\sigma}} \f$ is bounded by the three Rankine surfaces
   * \f$ f_A = \bar\sigma_A - f_\mathrm{t} \f$ of RankineSofteningPlasticity without hardening or softening
   * (associated flow, perfect plasticity, \f$ \Delta\boldsymbol{\varepsilon}^\mathrm{p} = \sum_A \Delta\lambda_A
   * \mathbf{E}_A \f$). Softening results from the scalar damage only,
   * \f[
   *   \boldsymbol{\sigma} = (1-\omega)\,\bar{\boldsymbol{\sigma}}, \qquad
   *   \omega = 1 - \exp(-\kappa_\mathrm{d}/\varepsilon_\mathrm{f}), \qquad
   *   \Delta\kappa_\mathrm{d} = \Delta\varepsilon^\mathrm{p}_\mathrm{V}/x_\mathrm{s}
   *     = \mathrm{tr}\,\Delta\boldsymbol{\varepsilon}^\mathrm{p} = \sum_A \Delta\lambda_A,
   * \f]
   * since all principal plastic strain increments are nonnegative for the Rankine surfaces (ductility measure
   * \f$ x_\mathrm{s} = 1 \f$). Hence \f$ \kappa_\mathrm{d} \f$ and \f$ \omega \f$ never decrease, and under uniaxial
   * tension \f$ \sigma = f_\mathrm{t} \exp(-\varepsilon^\mathrm{p}/\varepsilon_\mathrm{f}) \f$. The free energy is
   * \f$ \psi = (1-\omega)\,\bar\psi_\mathrm{el} \f$ with \f$ \bar\psi_\mathrm{el} = \frac12\bar{\boldsymbol{\sigma}}
   * :\mathbb{C}^{-1}:\bar{\boldsymbol{\sigma}} \f$, and the dissipation increment is \f$ (1-\omega)\,
   * \bar{\boldsymbol{\sigma}}:\Delta\boldsymbol{\varepsilon}^\mathrm{p} + \bar\psi_\mathrm{el}\,\Delta\omega \f$
   * (backward Euler). The consistent tangent reads \f$ (1-\omega)\,\bar{\mathbb{C}}^\mathrm{ep} -
   * \bar{\boldsymbol{\sigma}} \otimes \partial\omega/\partial\Delta\boldsymbol{\varepsilon} \f$ with
   * \f$ \partial\kappa_\mathrm{d}/\partial\Delta\boldsymbol{\varepsilon} = \boldsymbol{1} -
   * \bar{\mathbb{C}}^{\mathrm{ep},T}:\mathbb{C}^{-1}:\boldsymbol{1} \f$.
   *
   * Material properties: \f$ E, \nu, f_\mathrm{t}, \varepsilon_\mathrm{f} \f$ and optionally the density.
   * \f$ \varepsilon_\mathrm{f} \f$ is a material parameter (no regularization by the element size);
   * \f$ \varepsilon_\mathrm{f} = \infty \f$ recovers perfect Rankine plasticity.
   */
  class RankineDamagePlasticity : public MaterialHypoElastic {

  public:
    RankineDamagePlasticity( const double* materialProperties, const int nMaterialProperties, const int materialLabel );

    RankineDamagePlasticity( const RankineDamagePlasticity& )            = delete;
    RankineDamagePlasticity& operator=( const RankineDamagePlasticity& ) = delete;

    void computeStress( state3D&              state,
                        Math::Matrix6d&       dStressDDStrain,
                        const Math::Vector6d& dStrain,
                        const timeInfo&       timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    /// the stress update of computeStress, without the check of its results
    void updateStress( state3D&              state,
                       Math::Matrix6d&       dStressDDStrain,
                       const Math::Vector6d& dStrain,
                       const timeInfo&       timeInfo ) const;

    /// properties of the perfectly plastic RankineSofteningPlasticity for the effective stress
    const std::array< double, 6 >    effectiveProperties;
    const RankineSofteningPlasticity effectiveMaterial;
  };

} // namespace Marmot::Materials
