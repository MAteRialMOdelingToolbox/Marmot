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
#include "Marmot/MarmotPortability.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @brief Optional capability of a MaterialHypoElastic: a stress update that also yields the energy tangents.
   * @details A material provides the derivatives of its updated elastic energy density and dissipation wrt the strain
   * increment by inheriting from this class next to MaterialHypoElastic. Elements find the capability with
   * `dynamic_cast< const HypoElasticEnergyTangentProvider* >( material )` and call this variant instead of
   * MaterialHypoElastic::computeStress only when they assemble energy gradients, e.g., for dissipation based arc length
   * control. The standard material interface is not affected.
   */
  class MARMOT_API HypoElasticEnergyTangentProvider {

  public:
    virtual ~HypoElasticEnergyTangentProvider() = default;

    /**
     * @brief MaterialHypoElastic::computeStress, which also yields the energy tangents of the update.
     * @param[in,out] state Stress, energies and state variables, as for MaterialHypoElastic::computeStress.
     * @param[out] dStress_dStrain Algorithmic tangent.
     * @param[in] dStrain Linearized strain increment.
     * @param[in] timeInfo Time and time increment.
     * @param[out] dElasticEnergyDensity_dStrain Derivative of the updated elastic energy density wrt the strain
     * increment (Voigt, engineering shear strains).
     * @param[out] dDissipation_dStrain Derivative of the updated dissipation wrt the strain increment.
     */
    virtual void computeStressWithEnergyTangents( MaterialHypoElastic::state3D&        state,
                                                  Marmot::Math::Matrix6d&              dStress_dStrain,
                                                  const Marmot::Math::Vector6d&        dStrain,
                                                  const MaterialHypoElastic::timeInfo& timeInfo,
                                                  Marmot::Math::Vector6d&              dElasticEnergyDensity_dStrain,
                                                  Marmot::Math::Vector6d&              dDissipation_dStrain ) const = 0;
  };

} // namespace Marmot::Materials
