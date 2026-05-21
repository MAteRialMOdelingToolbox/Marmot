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
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include <Eigen/Core>
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Marmot::Materials {

  class VonMisesModel;

  class VonMisesInterface : public MarmotInterfaceMaterialHypoElastic {

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

    // Re-mapped properties for VonMisesModel: [E, nu, yieldStress, HLin, deltaYieldStress, delta]
    // Stored as a member so that VonMisesModel can hold a pointer to them for its lifetime
    std::array< double, 6 > vonMisesProps;

    // VonMisesModel is an implementation detail; keep it out of this public header.
    std::unique_ptr< VonMisesModel > vonMisesModel;

  public:
    VonMisesInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );
    ~VonMisesInterface() override;

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override;

    void initializeStateLayout() override;

    double getDensity() override;
  };

} // namespace Marmot::Materials
