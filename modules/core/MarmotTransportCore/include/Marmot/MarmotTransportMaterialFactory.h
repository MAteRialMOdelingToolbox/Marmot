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
#include "Marmot/MarmotPortability.h"
#include "Marmot/MarmotTransportMaterial.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace Marmot::Registry {

  /**
   * @class Marmot::Registry::TransportMaterialFactory
   * @brief Factory for transport materials, see Marmot::Materials::TransportMaterial.
   *
   * Transport materials register themselves by name, and transport elements create them by the name given in the
   * material section.
   */
  class MARMOT_API TransportMaterialFactory {
  public:
    /// Function signature for the factory functions registered in the map
    using transportMaterialFactoryFunction = std::function<
      Marmot::Materials::TransportMaterial*( const double* materialProperties,
                                             int           nMaterialProperties,
                                             int           materialNumber ) >;

    TransportMaterialFactory() = delete;

    /**
     * @brief Create a transport material by its name.
     * @param[in] name                Registered name of the material.
     * @param[in] materialProperties  Array of properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber      Identifier for the instance.
     * @return Pointer to the new material.
     * @throws std::invalid_argument if no transport material is registered under @p name.
     */
    static Marmot::Materials::TransportMaterial* createMaterial( const std::string& name,
                                                                 const double*      materialProperties,
                                                                 int                nMaterialProperties,
                                                                 int                materialNumber );

    /**
     * @brief Register a transport material with an auto-generated factory function.
     * @param[in] name Unique name under which the material is registered.
     * @return True if the registration was successful.
     */
    template < class T >
    static bool registerMaterial( const std::string& name )
    {
      auto& map = transportMaterialFactoryFunctionByName();

      assert( map.find( name ) == map.end() && "Transport material already registered!" );

      map[name] = []( const double* materialProperties,
                      int           nMaterialProperties,
                      int           materialNumber ) -> Marmot::Materials::TransportMaterial* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

  private:
    using TransportMaterialFactoryMap = std::unordered_map< std::string, transportMaterialFactoryFunction >;
    static TransportMaterialFactoryMap& transportMaterialFactoryFunctionByName();
  };
} // namespace Marmot::Registry
