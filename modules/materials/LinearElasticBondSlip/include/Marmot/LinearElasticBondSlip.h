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
#include "Marmot/MarmotBondSlipLaw.h"

namespace Marmot::Materials {

  /**
   * @class LinearElasticBondSlip
   * @brief Linear elastic bond-slip law.
   *
   * \f[ \tau_t = K_t\, s_t, \qquad \sigma_{n,i} = K_n\, s_{n,i} \f]
   *
   * With large stiffnesses it is a penalty formulation of perfect bond.
   *
   * Material properties:
   *  - @b K_t: tangential bond stiffness (stress per slip)
   *  - @b K_n: normal bond stiffness (stress per relative displacement)
   */
  class LinearElasticBondSlip : public Marmot::Materials::BondSlipLaw {
  public:
    /// tangential bond stiffness
    const double& Kt;
    /// normal bond stiffness
    const double& Kn;

    LinearElasticBondSlip( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeBondStress( State&                        state,
                            Marmot::Math::Matrix3d&       dBondStress_dSlip,
                            const Marmot::Math::Vector3d& slip,
                            const Marmot::Math::Vector3d& dSlip,
                            const TimeIncrement&          timeIncrement ) const override;
  };
} // namespace Marmot::Materials
