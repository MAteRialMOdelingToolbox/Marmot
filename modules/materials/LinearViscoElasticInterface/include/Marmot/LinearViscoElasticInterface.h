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
#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear elastic interface material
   * for 3D stress states.
   *
   * For further information see \ref linearelasticinterface.
   * according to the LinearViscoelasticPowerLaw model by Bazant et al. (2015)

   * generalized for 3D stress states.

   *

   * For further information see \ref b4.

   */
  class LinearViscoElasticInterface : public MarmotMaterialHypoElasticInterface {

    /// \brief Young's modulus
    const double& E_0;

    /// \brief Poisson's ratio
    const double& nu_0;

    /// \brief height of the middle layer
    const double& h;

    /// \brief power law compliance parameter for interphase layer
    const double& m;

    /// \brief power law exponent for interphase layer
    const double& n;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance for interphase layer
    const size_t nMaxwell;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain for interphase layer
    const double& minTau;

    /// \brief ratio of simulation time to days
    const double& timeToDays;

    class LinearViscoElasticInterfaceStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "MaxwellStateVars_force_uu", .length = 3 * 1 },
        { .name = "MaxwellStateVars_force_us", .length = 3 * 1 },
        { .name = "MaxwellStateVars_surface_stress_Z", .length = 9 * 1 },
        { .name = "MaxwellStateVars_surface_stress_Y", .length = 9 * 1 },
        { .name = "MaxwellStateVars_surface_stress_us", .length = 9 * 1 },
      } );

      WiechertInterface::mapStateVarMatrix_force_uu          MaxwellStateVars_force_uu;
      WiechertInterface::mapStateVarMatrix_force_us          MaxwellStateVars_force_us;
      WiechertInterface::mapStateVarMatrix_surface_stress_Z  MaxwellStateVars_surface_stress_Z;
      WiechertInterface::mapStateVarMatrix_surface_stress_Y  MaxwellStateVars_surface_stress_Y;
      WiechertInterface::mapStateVarMatrix_surface_stress_us MaxwellStateVars_surface_stress_us;
      LinearViscoElasticInterfaceStateVarManager( double* theStateVarVector, int nMaxwellUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          MaxwellStateVars_force_uu( &find( "MaxwellStateVars_force_uu" ), 3, nMaxwellUnits ),
          MaxwellStateVars_force_us( &find( "MaxwellStateVars_force_us" ), 3, nMaxwellUnits ),
          MaxwellStateVars_surface_stress_Z( &find( "MaxwellStateVars_surface_stress_Z" ), 9, nMaxwellUnits ),
          MaxwellStateVars_surface_stress_Y( &find( "MaxwellStateVars_surface_stress_Y" ), 9, nMaxwellUnits ),
          MaxwellStateVars_surface_stress_us( &find( "MaxwellStateVars_surface_stress_us" ), 9, nMaxwellUnits ){};
    };

    ::std::unique_ptr< LinearViscoElasticInterfaceStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface;
    using Tensor1D = Fastor::Tensor< double, 3 >;
    using Tensor2D = Fastor::Tensor< double, 3, 3 >;

    LinearViscoElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( double*       force,
                        double*       surfaceStress,
                        double*       H_inv_ij,
                        double*       Z_ijkl,
                        double*       H_inv_nF_ijk,
                        double*       Yn_H_inv_Fn_ijkl,
                        const double* dU,
                        const double* dSurfaceStrain,
                        const double* normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const ::std::string& stateName );

  private:
    WiechertInterface::Properties elasticModuli;
    WiechertInterface::Properties relaxationTimes;
    double                        zerothWiechertStiffness;

    static constexpr int powerLawApproximationOrder = 1;
  };
} // namespace Marmot::Materials
