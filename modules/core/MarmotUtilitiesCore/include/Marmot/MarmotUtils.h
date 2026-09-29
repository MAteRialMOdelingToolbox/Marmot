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
#include <stdexcept>
#include <string>

/** @struct StateView
 * @brief Structure to hold a pointer to the state location and its size.
 *
 * This structure is used to provide a view of the state variables in a finite element or material,
 * allowing access to the state data without copying it.
 */
struct StateView {
  double* stateLocation; ///< Pointer to the first element of the state variable block.
  int     stateSize;     ///< Number of `double` values in the state variable block.
};

namespace Marmot {

  /**
   * @brief Access a material property after checking that it has been provided.
   * @details Materials that bind their properties to references in the member initializer list cannot validate
   * the length of the property array in the constructor body first. Binding through this accessor turns a too short
   * property array into an error instead of a read beyond its end.
   * @param[in] materialProperties Array of the material properties.
   * @param[in] nMaterialProperties Length of @p materialProperties.
   * @param[in] index Index of the requested property.
   * @return Reference to the property.
   * @throws std::invalid_argument if @p index is not a valid index of @p materialProperties.
   */
  inline const double& checkedMaterialProperty( const double* materialProperties, int nMaterialProperties, int index )
  {
    if ( index < 0 || index >= nMaterialProperties )
      throw std::invalid_argument( "material property " + std::to_string( index ) + " is required, but only " +
                                   std::to_string( nMaterialProperties ) + " material properties are given" );
    return materialProperties[index];
  }

} // namespace Marmot
