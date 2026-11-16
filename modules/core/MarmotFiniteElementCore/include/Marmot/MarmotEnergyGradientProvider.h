/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
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
#include "Marmot/MarmotPortability.h"

namespace Marmot::Discretizations::Elements {

  /**
   * @brief Optional capability of an Element: the gradients of its elastic energy and dissipation wrt its dofs.
   * @details An element provides the gradients by inheriting from this class next to Element. Frameworks find the
   * capability with `dynamic_cast< EnergyGradientProvider* >( element )`, e.g., for dissipation based arc length
   * control, and assign the output arrays once. The element then writes the gradients in every computeKernels() call,
   * for the state of that call, consistent with its algorithmic update; computeKernelsExplicit() does not write them.
   * Elements without assigned outputs evaluate their kernels as before, and the standard element interface is not
   * affected.
   */
  class MARMOT_API EnergyGradientProvider {

  public:
    virtual ~EnergyGradientProvider() = default;

    /**
     * @brief Assign the arrays the element writes its energy gradients to in every computeKernels() call.
     * @param dElasticEnergy_dQ Gradient of the elastic energy, of size Element::getNDofPerElement(), in the element's
     * dof order; nullptr, together with @p dDissipation_dQ, stops the output.
     * @param dDissipation_dQ Gradient of the dissipation, of the same size and order.
     * @details Call it after the material section is assigned. It throws std::invalid_argument if the section or the
     * material cannot provide the gradients, so unsupported setups are rejected before the analysis. The arrays must
     * outlive the assignment.
     */
    virtual void assignEnergyGradientOutput( double* dElasticEnergy_dQ, double* dDissipation_dQ ) = 0;
  };

} // namespace Marmot::Discretizations::Elements
