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
#include "Marmot/MarmotBondSlipLaw.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace Marmot::Registry {

  /**
   * @class Marmot::Registry::BondSlipLawFactory
   * @brief Factory for bond-slip laws, see Marmot::Materials::BondSlipLaw.
   *
   * Bond-slip laws register themselves by name, and bond elements create them by the name given in the material
   * section.
   */
  class MARMOT_API BondSlipLawFactory {
  public:
    /// Function signature for the factory functions registered in the map
    using bondSlipLawFactoryFunction = std::function< Marmot::Materials::BondSlipLaw*( const double* materialProperties,
                                                                                       int nMaterialProperties,
                                                                                       int materialNumber ) >;

    BondSlipLawFactory() = delete;

    /**
     * @brief Create a bond-slip law by its name.
     * @param[in] name                Registered name of the bond-slip law.
     * @param[in] materialProperties  Array of properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber      Identifier for the instance.
     * @return Pointer to the new bond-slip law.
     * @throws std::invalid_argument if no bond-slip law is registered under @p name.
     */
    static Marmot::Materials::BondSlipLaw* createBondSlipLaw( const std::string& name,
                                                              const double*      materialProperties,
                                                              int                nMaterialProperties,
                                                              int                materialNumber );

    /**
     * @brief Register a bond-slip law with an auto-generated factory function.
     * @param[in] name Unique name under which the bond-slip law is registered.
     * @return True if the registration was successful.
     */
    template < class T >
    static bool registerBondSlipLaw( const std::string& name )
    {
      auto& map = bondSlipLawFactoryFunctionByName();

      assert( map.find( name ) == map.end() && "Bond-slip law already registered!" );

      map[name] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> Marmot::Materials::BondSlipLaw* { return new T( materialProperties, nMaterialProperties, materialNumber ); };
      return true;
    }

  private:
    using BondSlipLawFactoryMap = std::unordered_map< std::string, bondSlipLawFactoryFunction >;
    static BondSlipLawFactoryMap& bondSlipLawFactoryFunctionByName();
  };
} // namespace Marmot::Registry
