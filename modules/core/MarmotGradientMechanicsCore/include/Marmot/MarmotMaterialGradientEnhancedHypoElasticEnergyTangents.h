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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @brief Optional capability of a MaterialGradientEnhancedHypoElastic: a stress update that also yields the energy
   * tangents.
   * @tparam nNonlocalVariables Number of nonlocal variables of the material.
   * @details A material provides the derivatives of its updated elastic energy density and dissipation wrt the strain
   * increment and the nonlocal variables by inheriting from this class next to MaterialGradientEnhancedHypoElastic.
   * Elements find the capability with `dynamic_cast` and call this variant instead of
   * MaterialGradientEnhancedHypoElastic::computeStress only when they assemble energy gradients. The standard material
   * interface is not affected.
   */
  template < int nNonlocalVariables >
  class GradientEnhancedHypoElasticEnergyTangentProvider {

  public:
    using Material = MaterialGradientEnhancedHypoElastic< nNonlocalVariables >;

    /// @brief Derivatives of the updated elastic energy density and dissipation.
    struct energyTangents {
      /// @brief Elastic energy density wrt the strain increment (Voigt, engineering shear strains).
      Marmot::Math::Vector6d dElasticEnergyDensityddStrain = Marmot::Math::Vector6d::Zero();
      /// @brief Dissipation wrt the strain increment.
      Marmot::Math::Vector6d dDissipationddStrain = Marmot::Math::Vector6d::Zero();
      /// @brief Elastic energy density wrt the nonlocal variables.
      Eigen::Vector< double, nNonlocalVariables >
        dElasticEnergyDensityddK = Eigen::Vector< double, nNonlocalVariables >::Zero();
      /// @brief Dissipation wrt the nonlocal variables.
      Eigen::Vector< double, nNonlocalVariables > dDissipationddK = Eigen::Vector< double, nNonlocalVariables >::Zero();
    };

    virtual ~GradientEnhancedHypoElasticEnergyTangentProvider() = default;

    /**
     * @brief MaterialGradientEnhancedHypoElastic::computeStress, which also yields the energy tangents of the update.
     * @param[in,out] res Response, as for MaterialGradientEnhancedHypoElastic::computeStress.
     * @param[in,out] tan Algorithmic tangents.
     * @param[in] inc Increment of the strain and the nonlocal variables.
     * @param[out] energyTan Derivatives of the updated elastic energy density and dissipation.
     */
    virtual void computeStressWithEnergyTangents( typename Material::response&        res,
                                                  typename Material::tangents&        tan,
                                                  const typename Material::increment& inc,
                                                  energyTangents&                     energyTan ) const = 0;
  };

} // namespace Marmot::Materials
