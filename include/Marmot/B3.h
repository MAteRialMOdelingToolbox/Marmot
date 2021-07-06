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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/Solidification.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {

  class B3 : public MarmotMaterialHypoElastic {

    const double& nu;
    const double& q1;
    const double& q2;
    const double& q3;
    const double& q4;
    const double& q5;
    const double& hEnv;
    const double& shrinkageHalfTime;
    const double& ultimateShrinkageStrain;
    const double& n;
    const double& m;
    const size_t  numberOfKelvinUnits;
    const double& minimalRetardationTime;
    const size_t  numberOfKelvinUnitsDrying;
    const double& minimalRetardationTimeDrying;
    const double& dryingStart;
    const double& dTStatic;
    const double& timeToDays;
    const double& castTime;

    Solidification::MaterialParameters    basicCreepMaterialParameters;
    Solidification::KelvinChainProperties kelvinChainProperties;

    class B3StateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "EStatic", .length = 1 },
        { .name = "kelvinStateVars", .length = 0 },
      } );

      double&                               EStatic;
      Solidification::mKelvinStateVarMatrix kelvinStateVars;

      B3StateVarManager( double* theStateVarVector, int nKelvinUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          EStatic( find( "EStatic" ) ),
          kelvinStateVars( &find( "kelvinStateVars" ), 6, nKelvinUnits){};
    };
    std::unique_ptr< B3StateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    B3( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const std::string& stateName );
  };
} // namespace Marmot::Materials
