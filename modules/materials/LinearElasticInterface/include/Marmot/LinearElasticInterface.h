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
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear elastic interface material
   * for 3D stress states.
   *
   */
  class LinearElasticInterface : public MarmotInterfaceMaterialHypoElastic {
  public:
    using MarmotInterfaceMaterialHypoElastic::MarmotInterfaceMaterialHypoElastic;
    using Tensor1D = Marmot::FastorStandardTensors::Tensor3d;
    using Tensor2D = Marmot::FastorStandardTensors::Tensor33d;

    LinearElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void initializeStateLayout() override;

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override;

    int getNumberOfRequiredStateVars() const override { return 0; }

  private:
    /// \brief Young's modulus
    const double& E_0;

    /// \brief Poisson's ratio
    const double& nu_0;

    /// \brief height of the middle layer
    const double& h;
  };
} // namespace Marmot::Materials
