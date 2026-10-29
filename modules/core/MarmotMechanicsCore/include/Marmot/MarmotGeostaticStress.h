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

#include "Marmot/MarmotJournal.h"
#include <Eigen/Core>
#include <functional>

/**
 * @file MarmotGeostaticStress.h
 * @brief Utilities for evaluating geostatic (gravity-induced) initial stress states.
 *
 * Provides functions for computing the in-situ stress components from a
 * linearly-varying stress distribution defined along a vertical coordinate axis.
 */

namespace Marmot::ContinuumMechanics::GeostaticStress {

  /**
   * @brief Computes the three principal geostatic stress components at a given depth.
   *
   * The stress distribution is assumed to vary linearly with the vertical coordinate
   * @p coordinate_y as specified by @p geostaticStressDefintion.
   *
   * @param geostaticStressDefintion  Pointer to the array defining the linear distribution.
   * @param coordinate_y              Vertical coordinate of the evaluation point.
   * @return Tuple of three stress components (vertical, horizontal-1, horizontal-2).
   */
  std::tuple< double, double, double > getGeostaticStressFromLinearDistribution( const double* geostaticStressDefintion,
                                                                                 double        coordinate_y );

} // namespace Marmot::ContinuumMechanics::GeostaticStress
