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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotWiechert.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear visco elastic material
   * for 3D stress states.
   *
   * For further information see \ref linearviscoelasticwiechert.
   * according to the LinearViscoelasticPowerLaw model by Bazant et al. (2015)

   * generalized for 3D stress states.

   *

   * For further information see \ref b4.

   */
  class LinearViscoElasticWiechert : public MarmotMaterialHypoElastic {

    /// \brief Young's modulus
    const double& E;

    /// \brief Poisson's ratio
    const double& nu;

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

    class LinearViscoElasticWiechertStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "MaxwellStateVars", .length = 6 * 1 },
      } );

      Wiechert::mapStateVarMatrix MaxwellStateVars;
      LinearViscoElasticWiechertStateVarManager( double* theStateVarVector, int nMaxwellUnits )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          MaxwellStateVars( &find( "MaxwellStateVars" ), 6, nMaxwellUnits ){};
    };

    ::std::unique_ptr< LinearViscoElasticWiechertStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;
    using Tensor1D = Fastor::Tensor< double, 3 >;
    using Tensor2D = Fastor::Tensor< double, 3, 3 >;

    LinearViscoElasticWiechert( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void initializeStateLayout() override;

    void computeStress( state3D& state, double* C, const double* dStrain, const timeInfo& timeInfo ) const override;

    StateView getStateView( const ::std::string& stateName );

  private:
    Wiechert::Properties elasticModuli;
    Wiechert::Properties relaxationTimes;
    double               zerothWiechertStiffness;

    static constexpr int powerLawApproximationOrder = 1;

  private:
  };

} // namespace Marmot::Materials
