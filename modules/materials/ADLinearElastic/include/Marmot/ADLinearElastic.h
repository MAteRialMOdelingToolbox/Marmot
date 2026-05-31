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
#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot;

namespace Marmot::Materials {
  /**
   * @brief Implementation of a isotropic linear elastic material
   * for 3D stress states using automatic differentiation.
   *
   */
  class ADLinearElastic : public MarmotMaterialHypoElasticAD {
  public:
    using MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD;

    /// @brief Elasticity tensor for isotropic materials
    const autodiff::MatrixXdual C;

    ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber );

    double getDensity( const double* stateVars ) const override;

  protected:
    void computeStressAD( state3DAD&                 state,
                          const Marmot::Vector6dual& dStrain,
                          const timeInfo&            timeInfo ) const override;
  };
} // namespace Marmot::Materials
