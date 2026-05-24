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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @brief Factory class for creating MarmotMaterialGeneralGradientEnhancedHypoElastic instances.
   * This class uses the Factory Design Pattern to manage the creation of different
   * MarmotMaterialGeneralGradientEnhancedHypoElastic materials based on unique material names.
   */
  template < int nNonlocalVariables >
  class MarmotMaterialGeneralGradientEnhancedHypoElasticFactory {
  public:
    using materialFactoryFunction = std::function< MarmotMaterialGeneralGradientEnhancedHypoElastic<
      nNonlocalVariables >*( const double* materialProperties, int nMaterialProperties, int materialNumber ) >;

    MarmotMaterialGeneralGradientEnhancedHypoElasticFactory() = delete;

    /**
     * @brief Create a material instance based on its unique name and properties.
     * @param[in] materialName Unique name of the material to create.
     * @param[in] materialProperties Array of material properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber Unique identifier for the material instance.
     * @return Pointer to the created MarmotMaterial instance, or nullptr if creation failed.
     */
    static MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >* createMaterial(
      const std::string& materialName,
      const double*      materialProperties,
      int                nMaterialProperties,
      int                materialNumber );

    /**
     * @brief Register a material with its code and factory function.
     * @param[in] materialName Name of the material.
     * @return True if registration was successful, false if the code already exists.
     */
    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    /// @brief Get the map of material factory functions by material name.
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };
} // namespace MarmotLibrary
