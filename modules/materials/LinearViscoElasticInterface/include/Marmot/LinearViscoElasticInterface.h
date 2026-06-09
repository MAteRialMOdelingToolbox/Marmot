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
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear visco elastic interface material
   * for 3D stress states.
   *
   * For further information see \ref linearviscoelasticinterface.
   * according to the Wiechert model

   * generalized for 3D stress states.

   *

   * For further information see \ref b4.

   */
  class LinearViscoElasticInterface : public MarmotInterfaceMaterialHypoElastic {

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

  public:
    using MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic;

    LinearViscoElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override;

    int getNumberOfRequiredStateVars() const override;

  private:
    /// @brief Young's modulus of the #nKelvin Kelvin units
    WiechertInterface::Properties elasticModuli;
    /// @brief relaxation times of the #nKelvin Kelvin units
    WiechertInterface::Properties relaxationTimes;
    /// @brief stiffness of the zeroth Kelvin unit (the elastic response)
    double zerothWiechertStiffness;

    static constexpr int powerLawApproximationOrder = 1;
  };
} // namespace Marmot::Materials
