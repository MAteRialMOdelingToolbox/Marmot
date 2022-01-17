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
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotSolidification.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {

  class B4 : public MarmotMaterialHypoElastic {

    // elasticity
    const double& nu;

    // basic creep
    const double& q1;
    const double& q2;
    const double& q3;
    const double& q4;
    const double& n;
    const double& m;
    const size_t  nKelvinBasic;
    const double& minTauBasic;

    // autogenous shrinkage
    const double& ultimateAutogenousShrinkageStrain;
    const double& autogenousShrinkageHalfTime;
    const double& alpha;
    const double& rt;

    // drying shrinkage
    const double& ultimateDryingShrinkageStrain;
    const double& dryingShrinkageHalfTime;
    const double& dryingStart;
    const double& hEnv;

    // drying creep
    const double& q5;
    const size_t  nKelvinDrying;
    const double& minTauDrying;

    // time parameters
    const double& castTime;
    const double& timeToDays;

    class B4StateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kelvinStateVars", .length = 0 },
      } );

      KelvinChain::mapStateVarMatrix kelvinStateVars;

      B4StateVarManager( double* theStateVarVector, int nKelvinUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          kelvinStateVars( &find( "kelvinStateVars" ), 6, nKelvinUnits ){};
    };
    std::unique_ptr< B4StateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    B4( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const std::string& stateName );

  private:
    SolidificationTheory::Parameters            solidificationParameters;
    SolidificationTheory::KelvinChainProperties solidificationKelvinProperties;

    KelvinChain::Properties basicCreepElasticModuli;
    KelvinChain::Properties basicCreepRetardationTimes;

    static constexpr int dryingCreepComplianceApproximationOrder = 5;
    static constexpr int basicCreepComplianceApproximationOrder  = 2;

    template < typename T_ >
    T_ phi( T_ xi, double b, double xiZero )
    {
      T_ val = sqrt( exp( tanh( sqrt( xi - xiZero ) ) * b ) - exp( tanh( sqrt( -xiZero ) ) * b ) );
      return val;
    }
  };
} // namespace Marmot::Materials
