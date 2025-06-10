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
#include "autodiff/forward/dual.hpp"
#include <Eigen/src/Core/Map.h>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot;

namespace Marmot::Materials {
  /**
   * \brief Implementation of a isotropic linear elastic material
   * for 3D stress states using automatic differentiation.
   *
   */
  class ADLinearElastic : public MarmotMaterialHypoElasticAD {
  public:
    using MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD;

    const double& E;
    const double& nu;

    ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    void computeStressAD( autodiff::dual*       stress,
                          const autodiff::dual* dStrain,
                          const double*         timeOld,
                          const double          dT,
                          double&               pNewDT );

    StateView getStateView( const std::string& result ) { return { nullptr, 0 }; };

    int getNumberOfRequiredStateVars() { return 0; }
  };
} // namespace Marmot::Materials
