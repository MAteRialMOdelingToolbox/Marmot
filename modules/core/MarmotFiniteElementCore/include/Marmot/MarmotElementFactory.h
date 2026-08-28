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
#include "Marmot/MarmotElement.h"
#include "Marmot/marmot_export.h"
#include <cassert>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @class MarmotElementFactory
   * @brief Factory class for creating element instances.
   * This class provides a mechanism to register elements by their code and name,
   * and to create element instances based on their properties.
   */
  class MARMOT_EXPORT MarmotElementFactory {
  public:
    /// @brief Factory function pointer type: takes an element number and returns a new MarmotElement.
    using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
    MarmotElementFactory()       = delete;

    /**
     * @brief Create an element instance based on its name and number.
     * @param[in] elementName   Registered name of the element type.
     * @param[in] elementNumber Unique identifier for the element instance.
     * @return Pointer to the created MarmotElement instance.
     * @throws std::invalid_argument if @p elementName is not registered.
     */
    static MarmotElement* createElement( const std::string& elementName, int elementNumber )
    {
      auto& map = elementFactoryFunctionByName();
      auto  it  = map.find( elementName );
      if ( it == map.end() ) {
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << "Element " + elementName + " not registered!" );
      }

      return it->second( elementNumber );
    }

    /**
     * @brief Register an element with its name.
     * @param[in] elementName     Name of the element type to register.
     * @param[in] factoryFunction Factory function pointer that creates the element.
     * @return True if registration was successful.
     */
    static bool registerElement( const std::string& elementName, elementFactoryFunction factoryFunction )
    {
      auto& map = elementFactoryFunctionByName();

      assert( map.find( elementName ) == map.end() && "Element already registered!" );

      map[elementName] = factoryFunction;
      return true;
    }

  private:
    using ElementFactoryMap = std::unordered_map< std::string, elementFactoryFunction >;
    static ElementFactoryMap& elementFactoryFunctionByName();
  };

} // namespace MarmotLibrary
