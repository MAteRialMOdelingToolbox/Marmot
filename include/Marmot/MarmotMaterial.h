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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotUtils.h"
#include <string>

/**
 * @class MarmotMaterial
 * @brief Abstract base class for material models in the Marmot framework.
 *
 * Provides a common interface for material models, including access to material
 * properties, state variables, initialization, and density. Concrete material
 * implementations must override the pure virtual functions.
 */
class MarmotMaterial {

protected:
  const double* materialProperties; ///< Pointer to array of material property values.
  const int     nMaterialProperties; ///< Number of material properties.

  double* stateVars; ///< Pointer to array of state variables.
  int     nStateVars; ///< Number of assigned state variables.

public:
  const int materialNumber; ///< Identifier for material type/implementation.

  /**
   * @brief Construct a material object.
   * @param materialProperties Pointer to array of material properties.
   * @param nMaterialProperties Number of material properties.
   * @param materialNumber Unique identifier for the material.
   */
  MarmotMaterial( const double* materialProperties, int nMaterialProperties, int materialNumber );

  /** @brief Virtual destructor for safe polymorphic cleanup. */
  virtual ~MarmotMaterial();

  /**
   * @return Number of state variables required by the material.
   */
  virtual int getNumberOfRequiredStateVars() = 0;

  /**
   * @brief Assign state variable array to material.
   * @param stateVars Pointer to state variable array.
   * @param nStateVars Number of state variables.
   */
  virtual void assignStateVars( double* stateVars, int nStateVars );

  /**
   * @brief Access material state variables by name.
   * @param stateName Name of the requested state variable.
   * @return A view into the state variable array.
   */
  virtual StateView getStateView( const std::string& stateName ) = 0;

  /**
   * @brief Get pointer to currently assigned state variables.
   * @return Pointer to state variable array.
   */
  double* getAssignedStateVars();

  /**
   * @brief Get number of assigned state variables.
   * @return Number of assigned state variables.
   */
  int getNumberOfAssignedStateVars();

  /**
   * @brief Initialize material state (default: no action).
   */
  virtual void initializeYourself();

  /**
   * @brief Get material density.
   * @return Density value (default: 0 if not overridden).
   */
  virtual double getDensity();
};
