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

  /**
   * \brief Implementation of a linear elastic material
   * according to the LinearViscoelasticPowerLaw model by Bazant et al. (2015)
   * generalized for 3D stress states.
   *
   * For further information see \ref b4.
   */
  class LinearViscoelasticPowerLaw : public MarmotMaterialHypoElastic {

    /// \brief Young's modulus
    const double& E;

    /// \brief Poisson's ratio 
    const double& nu;

    /// \brief power law compliance parameter
    const double& m;

    /// \brief power law exponent
    const double& n;
    
    /// \brief number of Kelvin units to approximate the viscoelastic compliance
    const size_t nKelvin;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain
    const double& minTau;

    /// \brief ratio of simulation time to days
    const double& timeToDays;
    
    class LinearViscoelasticPowerLawStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kelvinStateVars", .length = 0 },
      } );

      KelvinChain::mapStateVarMatrix kelvinStateVars;

      LinearViscoelasticPowerLawStateVarManager( double* theStateVarVector, int nKelvinUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          kelvinStateVars( &find( "kelvinStateVars" ), 6, nKelvinUnits ){};
    };
    std::unique_ptr< LinearViscoelasticPowerLawStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearViscoelasticPowerLaw( const double* materialProperties, int nMaterialProperties, int materialLabel );

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

    KelvinChain::Properties elasticModuli;
    KelvinChain::Properties retardationTimes;

    static constexpr int powerLawApproximationOrder  = 5;

    /// \brief creep compliance function
    template < typename T_ >
    T_ phi( T_ tau, double m, double n )
    {
      T_ val = m * pow( tau, n );
      return val;
    }
  };
} // namespace Marmot::Materials
