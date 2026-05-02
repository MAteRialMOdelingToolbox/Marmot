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
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/VonMises.h"
#include "Marmot/VonMisesConstants.h"
#include <Eigen/Core>
#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {

  class VonMisesInterface : public MarmotMaterialHypoElasticInterface {

    /// \brief Young's modulus
    const double& E_0;

    /// \brief Poisson's ratio
    const double& nu_0;

    /// \brief height of the middle layer
    const double& h;

    // plasticity parameters
    const double& yieldStress;
    const double& HLin;
    const double& deltaYieldStress;
    const double& delta;

    class VonMisesInterfaceStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kappa", .length = 1 },
        { .name = "C_ep_voigt", .length = 36 },
      } );

      double&                                                      kappa;
      Eigen::Map< Eigen::Matrix< double, 6, 6, Eigen::RowMajor > > C_ep_voigt;

      VonMisesInterfaceStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          kappa( find( "kappa" ) ),
          C_ep_voigt( &find( "C_ep_voigt" ) ){};
    };
    std::unique_ptr< VonMisesInterfaceStateVarManager > managedStateVars;

    // Re-mapped properties for VonMisesModel: [E, nu, yieldStress, HLin, deltaYieldStress, delta]
    // Stored as a member so that VonMisesModel can hold a pointer to them for its lifetime
    std::array< double, 6 > vonMisesProps;

    // VonMisesModel instance — instantiated once in the constructor and reused every increment
    VonMisesModel vonMisesModel;

  public:
    VonMisesInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( double*       force,
                        double*       surfaceStress,
                        double*       Q_ij,
                        double*       Z_ijkl,
                        double*       H_ijk,
                        double*       Y_ijkl,
                        const double* dU,
                        const double* dSurfaceStrain,
                        const double* normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT ) override;

    int getNumberOfRequiredStateVars() override { return VonMisesInterfaceStateVarManager::layout.nRequiredStateVars; }

    void assignStateVars( double* stateVars, int nStateVars ) override;

    StateView getStateView( const std::string& stateName ) override;

    double getDensity() override;
  };

} // namespace Marmot::Materials
