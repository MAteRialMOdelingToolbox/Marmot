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
#include "Marmot/MarmotCell.h"

/**
 * @class MarmotCellElement
 * @brief A MarmotCell that defines its own material points, like a finite element defines its quadrature points.
 *
 * @details In contrast to a plain MarmotCell, which integrates over whatever material points currently lie
 * inside it, a cell element requests a fixed set of material points: the host framework queries
 * getNMaterialPoints(), getRequestedMaterialPointCoordinates() and getRequestedMaterialPointVolumes(), creates
 * one material point per request (e.g. at the points of the quadrature rule passed to
 * MarmotLibrary::MarmotCellElementFactory), and assigns them back with MarmotCell::assignMaterialPoints().
 * The cell element may additionally carry its own state variables.
 *
 * Instances are created by name via MarmotLibrary::MarmotCellElementFactory. No concrete cell element is
 * contained in Marmot at present.
 */
class MarmotCellElement : public MarmotCell {

public:
  /// Virtual destructor; cell elements are owned through base class pointers.
  virtual ~MarmotCellElement() = default;

  /**
   * @brief Number of material points requested by the cell element.
   * @return Number of material points.
   */
  virtual int getNMaterialPoints() const = 0;

  /**
   * @brief Coordinates of the requested material points.
   * @param[out] coordinates getNMaterialPoints() @f$ \times @f$ nDim values, point by point.
   */
  virtual void getRequestedMaterialPointCoordinates( double* coordinates ) const = 0;

  /**
   * @brief Volumes (integration weights) of the requested material points.
   * @param[out] volumes getNMaterialPoints() values.
   */
  virtual void getRequestedMaterialPointVolumes( double* volumes ) const = 0;

  /**
   * @brief Number of state variables required by the cell element itself.
   * @return Number of state variables.
   */
  virtual int getNumberOfRequiredStateVars() = 0;

  /**
   * @brief Assigns the (externally owned) state variable array of the cell element.
   * @param[in,out] stateVars  State variable array.
   * @param[in]     nStateVars Size of @p stateVars.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars ) = 0;
};
