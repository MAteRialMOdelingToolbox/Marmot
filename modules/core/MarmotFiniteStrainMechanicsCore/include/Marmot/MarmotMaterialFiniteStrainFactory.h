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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/marmot_export.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @class MarmotMaterialFiniteStrainFactory
   * @brief Factory class for creating finite-strain material instances.
   *
   * This class provides a mechanism to register materials by their name,
   * and to create material instances based on their properties.
   * It allows for dynamic material creation without hardcoding specific material types.
   */
  class MARMOT_EXPORT MarmotMaterialFiniteStrainFactory {
  public:
    /// Factory function type: creates a MarmotMaterialFiniteStrain from properties and a material number.
    using materialFactoryFunction = std::function<
      MarmotMaterialFiniteStrain*( const double* materialProperties, int nMaterialProperties, int materialNumber ) >;

    MarmotMaterialFiniteStrainFactory() = delete;

    /**
     * @brief Create a material instance based on its name and properties.
     * @param[in] materialName       Registered name of the material type.
     * @param[in] materialProperties Array of material properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber Unique identifier for the material instance.
     * @return Pointer to the created MarmotMaterialFiniteStrain instance.
     */
    static MarmotMaterialFiniteStrain* createMaterial( const std::string& materialName,
                                                       const double*      materialProperties,
                                                       int                nMaterialProperties,
                                                       int                materialNumber );

    /**
     * @brief Register a material type with its name.
     * @tparam T  Concrete material class to register; must be constructible from
     *            (const double* materialProperties, int nMaterialProperties, int materialNumber).
     * @param[in] materialName Registered name for the material type.
     * @return True if registration was successful.
     */
    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> MarmotMaterialFiniteStrain* { return new T( materialProperties, nMaterialProperties, materialNumber ); };
      return true;
    }

  private:
    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };
} // namespace MarmotLibrary
