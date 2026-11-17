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
    /// @brief The material base class, from which the class providing the capability derives as well.
    using MaterialBase = MaterialGradientEnhancedHypoElastic< nNonlocalVariables >;

    /// @brief Derivatives of the updated elastic energy density and dissipation.
    struct energyTangents {
      /// @brief Elastic energy density wrt the strain increment (Voigt, engineering shear strains).
      Marmot::Math::Vector6d dElasticEnergyDensity_dStrain = Marmot::Math::Vector6d::Zero();
      /// @brief Dissipation wrt the strain increment.
      Marmot::Math::Vector6d dDissipation_dStrain = Marmot::Math::Vector6d::Zero();
      /// @brief Elastic energy density wrt the nonlocal variables.
      Eigen::Vector< double, nNonlocalVariables >
        dElasticEnergyDensity_dK = Eigen::Vector< double, nNonlocalVariables >::Zero();
      /// @brief Dissipation wrt the nonlocal variables.
      Eigen::Vector< double, nNonlocalVariables > dDissipation_dK = Eigen::Vector< double, nNonlocalVariables >::Zero();
    };

    virtual ~GradientEnhancedHypoElasticEnergyTangentProvider() = default;

    /**
     * @brief MaterialGradientEnhancedHypoElastic::computeStress, which also yields the energy tangents of the update.
     * @param[in,out] res Response, as for MaterialGradientEnhancedHypoElastic::computeStress.
     * @param[in,out] tan Algorithmic tangents.
     * @param[in] inc Increment of the strain and the nonlocal variables.
     * @param[out] energyTan Derivatives of the updated elastic energy density and dissipation.
     */
    virtual void computeStressWithEnergyTangents( typename MaterialBase::response&        res,
                                                  typename MaterialBase::tangents&        tan,
                                                  const typename MaterialBase::increment& inc,
                                                  energyTangents&                         energyTan ) const = 0;

    /**
     * @brief MaterialGradientEnhancedHypoElastic::computePlaneStress, which also yields the energy tangents of the
     * converged 3D update.
     * @param[in,out] res Response, as for MaterialGradientEnhancedHypoElastic::computePlaneStress.
     * @param[in,out] tan Algorithmic tangents of the converged 3D update.
     * @param[in] inc Increment of the in-plane strain and the nonlocal variables.
     * @param[out] energyTan Derivatives of the updated elastic energy density and dissipation wrt the 3D strain
     * increment and the nonlocal variables; the condensation of the out-of-plane strain is left to the element.
     * @details The default is the plane stress iteration of MaterialGradientEnhancedHypoElastic::computePlaneStress
     * (MaterialGradientEnhancedHypoElastic::iteratePlaneStress) around computeStressWithEnergyTangents(). A material
     * that overrides MaterialGradientEnhancedHypoElastic::computePlaneStress overrides this method consistently, so
     * that the stress does not depend on whether energy tangents are requested. The class providing the capability
     * must also be a MaterialGradientEnhancedHypoElastic.
     */
    virtual void computePlaneStressWithEnergyTangents( typename MaterialBase::response&        res,
                                                       typename MaterialBase::tangents&        tan,
                                                       const typename MaterialBase::increment& inc,
                                                       energyTangents&                         energyTan ) const
    {
      const auto& material = dynamic_cast< const MaterialBase& >( *this );
      MaterialBase::iteratePlaneStress( material.getNumberOfRequiredStateVars(),
                                        res,
                                        tan,
                                        inc,
                                        StressUpdateWithEnergyTangents{ *this, energyTan } );
    }

  private:
    /// @brief computeStressWithEnergyTangents() as the 3D stress update of MaterialBase::iteratePlaneStress.
    struct StressUpdateWithEnergyTangents {
      /// @brief The material.
      const GradientEnhancedHypoElasticEnergyTangentProvider& provider;
      /// @brief The energy tangents of the last update, i.e., of the converged 3D state.
      energyTangents& energyTan;

      /**
       * @brief The 3D stress update with the energy tangents.
       * @param[in,out] res Response.
       * @param[in,out] tan Algorithmic tangents.
       * @param[in] inc Increment.
       */
      void operator()( typename MaterialBase::response&        res,
                       typename MaterialBase::tangents&        tan,
                       const typename MaterialBase::increment& inc ) const
      {
        provider.computeStressWithEnergyTangents( res, tan, inc, energyTan );
      }
    };
  };

} // namespace Marmot::Materials
