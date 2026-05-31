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
#include "Marmot/MarmotTypedefs.h"

/**
 * @file MarmotUtility.h
 * @brief Miscellaneous utility helpers for controlling increment sizes.
 */

namespace Marmot {
  /**
   * @brief Signals that the current increment should be discarded.
   *
   * Sets @p pNewDT to the requested @p value and prints a warning message.
   * Intended to be called inside material routines when convergence cannot
   * be achieved, instructing the host code to restart with a smaller step.
   *
   * @param pNewDT   Output parameter: suggested new time-step size.
   * @param value    Suggested value for the new time-step size.
   * @param message  Human-readable reason for discarding the increment.
   */
  void discardTheIncrement( double& pNewDT, double value, const std::string& message );
} // namespace Marmot
