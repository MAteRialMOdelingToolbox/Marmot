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
#include "Marmot/marmot_export.h"
#include <string>

/** @struct MarmotMaterialSection
 * @brief Structure to hold material section properties.
 *
 * This structure is used to define a material section with its code and properties,
 * allowing for flexible material definitions in finite element analysis.
 */
class MARMOT_EXPORT MarmotMaterialSection {
public:
  const std::string materialName;        ///< Name identifying the material model.
  const double*     materialProperties;  ///< Pointer to the array of material property values.
  int               nMaterialProperties; ///< Number of material property values.

  /**
   * @brief Construct a MarmotMaterialSection.
   * @param[in] materialName       Name identifying the material model.
   * @param[in] materialProperties Pointer to the array of material property values.
   * @param[in] nMaterialProperties Number of material property values.
   */
  MarmotMaterialSection( const std::string materialName, const double* materialProperties, int nMaterialProperties )
    : materialName( materialName ),
      materialProperties( materialProperties ),
      nMaterialProperties( nMaterialProperties ){};
};

/** @struct ElementProperties
 * @brief Structure to hold element properties.
 *
 * This structure is used to define properties of a finite element,
 * allowing for flexible element definitions in finite element analysis.
 */
class MARMOT_EXPORT ElementProperties {
public:
  const double* elementProperties;  ///< Pointer to the array of element property values.
  int           nElementProperties; ///< Number of element property values.

  /**
   * @brief Construct an ElementProperties object.
   * @param[in] elementProperties  Pointer to the array of element property values.
   * @param[in] nElementProperties Number of element property values.
   */
  ElementProperties( const double* elementProperties, int nElementProperties )
    : elementProperties( elementProperties ), nElementProperties( nElementProperties ){};
};
