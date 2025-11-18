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
#include "Marmot/MarmotElement.h"
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
  class MarmotElementFactory {
  public:
    using elementFactoryFunction = MarmotElement* (*)( int elementNumber );
    MarmotElementFactory()       = delete;

    /**
     * @brief Create an element instance based on its code and number.
     * @param[in] elementCode Unique code for the element.
     * @param[in] elementNumber Unique identifier for the element instance.
     * @return Pointer to the created MarmotElement instance, or nullptr if creation failed.
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
     * @param[in] elementName Name of the element.
     * @return True if registration was successful, false if the code already exists.
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
