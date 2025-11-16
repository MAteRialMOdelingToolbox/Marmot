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
  const double* materialProperties;  ///< Pointer to array of material property values.
  const int     nMaterialProperties; ///< Number of material properties.

  double* stateVars  = nullptr;      ///< Pointer to array of state variable values.
  int     nStateVars = 0;            ///< Number of state variables.

public:
  const int materialNumber; ///< Identifier for material type/implementation.

  /**
   * @brief Construct a material object.
   * @param[in] materialProperties Pointer to array of material properties.
   * @param[in] nMaterialProperties Number of material properties.
   * @param[in] materialNumber Unique identifier for the material.
   */
  MarmotMaterial( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : materialProperties( materialProperties ),
      nMaterialProperties( nMaterialProperties ),
      materialNumber( materialNumber )
  {
  }

  /** @brief Virtual destructor for safe polymorphic cleanup. */
  virtual ~MarmotMaterial() = default;
};
