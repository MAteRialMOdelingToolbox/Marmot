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
 * Thomas Mader  thomas.mader@uibk.ac.at
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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
  namespace ContinuumMechanics::LocalizationAnalysis {

    Marmot::Matrix3d computeAcousticTensor( const Marmot::Matrix6d& materialTangent,
                                            const Marmot::Vector3d& normalVector );

    bool localizationChecker( const Marmot::Matrix3d& acousticTensor );

    Marmot::Vector3d computeNormalVector( double alpha, double beta );

    double minimumDeterminantAcousticTensor( const Marmot::Matrix6d& materialTangent );

  } // namespace ContinuumMechanics::LocalizationAnalysis
} // namespace Marmot
