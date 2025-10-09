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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @brief An implementation of classical J2 plasticity with isotropic hardening.
   */
  class VonMisesModel : public MarmotMaterialHypoElastic {

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    void computeStress( double* stress,
                        double* dStress_dStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT ) override;

    /**
     * @brief Get material density.
     * @return Density value.
     * @throw std::runtime_error if density is not defined.
     */
    double getDensity() override;

    class VonMisesModelStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kappa", .length = 1 },
      } );

      /**
       * @brief Hardening parameter.
       */
      double& kappa;

      VonMisesModelStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), kappa( find( "kappa" ) ){};
    };
    std::unique_ptr< VonMisesModelStateVarManager > managedStateVars;

    int getNumberOfRequiredStateVars() override { return VonMisesModelStateVarManager::layout.nRequiredStateVars; }

    void assignStateVars( double* stateVars, int nStateVars ) override;

    StateView getStateView( const std::string& result ) override;
  };

} // namespace Marmot::Materials
