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
#include "Marmot/MarmotTransportMaterialFactory.h"
#include "Marmot/MarmotJournal.h"
#include <stdexcept>

Marmot::Materials::TransportMaterial* Marmot::Registry::TransportMaterialFactory::createMaterial(
  const std::string& name,
  const double*      materialProperties,
  int                nMaterialProperties,
  int                materialNumber )
{
  auto& map = transportMaterialFactoryFunctionByName();
  auto  it  = map.find( name );
  if ( it == map.end() ) {
    std::string registered = "Registered transport materials are: ";
    for ( const auto& pair : map )
      registered += pair.first + ", ";
    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__
                                 << " Transport material " + name + " not registered! " + registered );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}

Marmot::Registry::TransportMaterialFactory::TransportMaterialFactoryMap& Marmot::Registry::TransportMaterialFactory::
  transportMaterialFactoryFunctionByName()
{
  static TransportMaterialFactoryMap map;
  return map;
}
