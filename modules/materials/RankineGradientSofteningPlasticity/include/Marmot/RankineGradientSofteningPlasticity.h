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
#include "Marmot/MarmotMaterialGradientEnhancedHypoElasticEnergyTangents.h"
#include "Marmot/RankineSoftening.h"

namespace Marmot::Materials {

  /**
   * @brief Over-nonlocal (implicit gradient) version of RankineSofteningPlasticity: multi-surface Rankine plasticity
   * whose yield stress depends on a combination of the local and the nonlocal softening variable, following Vermeer
   * and Brinkgreve (1994) and the structure of RankineGradientDamagePlasticity.
   *
   * @details The three Rankine surfaces \f$ f_A = \sigma_A - f_\mathrm{y}(\hat\kappa) \f$ with associated flow and
   * \f$ \Delta\kappa = \sum_A \Delta\lambda_A \f$ are those of RankineSofteningPlasticity. The local softening variable
   * \f$ \kappa \f$ is the source of the Helmholtz equation solved by the element,
   * \f[
   *   \bar\kappa - l^2\,\Delta\bar\kappa = \kappa, \qquad
   *   \nabla\bar\kappa\cdot\boldsymbol{n} = 0 \text{ on the boundary},
   * \f]
   * and the yield stress follows from the over-nonlocal combination (\f$ m = 1 \f$: standard nonlocal, \f$ m > 1 \f$:
   * over-nonlocal, \f$ m = 0 \f$: local, i.e., RankineSofteningPlasticity)
   * \f[
   *   \hat\kappa = m\,\bar\kappa + (1-m)\,\kappa, \qquad
   *   f_\mathrm{y}(\hat\kappa) = \max\left(f_\mathrm{y}^0 + H_\mathrm{iso,lin}\hat\kappa
   *     + \Delta f_\mathrm{y}^{0,\infty}\left(1-\exp(-\delta\hat\kappa)\right),\ 0\right).
   * \f]
   * For a softening law and \f$ m > 1 \f$, the yield stress increases with \f$ \kappa \f$ at a fixed nonlocal variable,
   * so that the return mapping hardens locally and the plastic strain cannot localize in a set of measure zero. The
   * yield stress is floored at zero in terms of \f$ \hat\kappa \f$; the stress power \f$ \boldsymbol{\sigma}:
   * \Delta\boldsymbol{\varepsilon}^\mathrm{p} = f_\mathrm{y}(\hat\kappa)\,\Delta\kappa \f$ is the dissipation
   * increment, the free energy is the elastic strain energy.
   *
   * The backward Euler return mapping in principal space is that of RankineSofteningPlasticity with the hardening
   * modulus \f$ (1-m)\,f_\mathrm{y}' \f$. Returned to the element: \f$ \kappa \f$ (KLocal), \f$ c = l^2 \f$, the
   * consistent tangent \f$ \partial\boldsymbol{\sigma}/\partial\Delta\boldsymbol{\varepsilon} \f$,
   * \f$ \partial\boldsymbol{\sigma}/\partial\bar\kappa = -\sum_A (\mathbb{C}_\mathrm{p}\,\partial\Delta\boldsymbol{
   * \lambda}/\partial\bar\kappa)_A\,\boldsymbol{n}_A\otimes\boldsymbol{n}_A \f$ with \f$ \partial\Delta\boldsymbol{
   * \lambda}/\partial\bar\kappa = -m\,f_\mathrm{y}'\,\boldsymbol{J}^{-1}\boldsymbol{1} \f$ on the active surfaces,
   * \f$ \partial\kappa/\partial\Delta\boldsymbol{\varepsilon} = \boldsymbol{1} - (\partial\boldsymbol{\sigma}/\partial
   * \Delta\boldsymbol{\varepsilon})^T:\mathbb{C}^{-1}:\boldsymbol{1} \f$ and \f$ \partial\kappa/\partial\bar\kappa =
   * \sum_A \partial\Delta\lambda_A/\partial\bar\kappa \f$. Since the yield stress depends on the nonlocal variable, a
   * step without a strain increment may be plastic.
   *
   * Material properties: \f$ E, \nu, f_\mathrm{y}^0, H_\mathrm{iso,lin}, \Delta f_\mathrm{y}^{0,\infty}, \delta, l, m
   * \f$, optionally the density and the nonlocal viscosity (default zero, for explicit dynamics only).
   */
  class RankineGradientSofteningPlasticity : public MaterialGradientEnhancedHypoElastic< 1 >,
                                             public GradientEnhancedHypoElasticEnergyTangentProvider< 1 > {

  public:
    using Base = MaterialGradientEnhancedHypoElastic< 1 >;
    using typename Base::increment;
    using typename Base::response;
    using typename Base::tangents;

    RankineGradientSofteningPlasticity( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    /**
     * @brief computeStress(), which also yields the derivatives of the updated elastic energy density and
     * dissipation wrt the strain increment and the nonlocal variable.
     * @param[in,out] res As for computeStress().
     * @param[out] tan As for computeStress().
     * @param[in] inc As for computeStress().
     * @param[out] energyTan The energy tangents.
     */
    void computeStressWithEnergyTangents( response&        res,
                                          tangents&        tan,
                                          const increment& inc,
                                          energyTangents&  energyTan ) const override;

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    /// the stress update with the check of its results; energy tangents only if energyTan is given
    void computeStressImplicit( response& res, tangents& tan, const increment& inc, energyTangents* energyTan ) const;

    /// the stress update, without the check of its results
    void updateStress( response& res, tangents& tan, const increment& inc, energyTangents* energyTan ) const;

    const RankineSoftening::Hardening hardening;
  };

} // namespace Marmot::Materials
