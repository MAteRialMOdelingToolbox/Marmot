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

#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"

#include <memory>
#include <vector>

namespace Marmot::Materials {

  class LinearViscoelasticPowerLaw;

  /**
   * @brief Interface-material wrapper around LinearViscoelasticPowerLaw.
   *
   * The interface kinematics are converted to an average three-dimensional
   * strain increment. The wrapped bulk material performs the Kelvin-chain
   * stress update, after which the stress and tangent are reduced to the
   * interface force, surface stress, and interface tangent operators.
   *
   * Material properties are ordered as
   * `[E, nu, h, m, n, nKelvin, minTau, timeToDays, optional density]`.
   */
  class KelvinChainInterface : public MarmotInterfaceMaterialHypoElastic {
  public:
    KelvinChainInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );
    ~KelvinChainInterface() override;

    void computeStress( State&               state,
                        Tangents&            tangents,
                        const Deformation&   deformation,
                        const TimeIncrement& timeIncrement ) override;

    void initializeYourself( double* stateVars, int nStateVars ) override;

    double getDensity() override;

  private:
    const double& h;

    // Properties passed to the wrapped bulk material, with interface thickness removed.
    std::vector< double >                         bulkMaterialProperties;
    std::unique_ptr< LinearViscoelasticPowerLaw > bulkMaterial;
  };

} // namespace Marmot::Materials
