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
#include "Marmot/MarmotUtils.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

/// @brief A convenience auxiliary class for managing multiple statevars with arbitrary length in a single consecutive
/// double array
class MarmotStateVarVectorManager {

public:
  /// get a StateView for a statevar entry
  inline StateView getStateView( const std::string& name ) const
  {
    auto entry = theLayout.entries.at( name );
    return { theStateVars + entry.index, entry.length };
  }

  /// get the reference to the first array element of an entry in the statevar vector
  inline double& find( const std::string& name ) const { return theStateVars[theLayout.entries.at( name ).index]; }

  /// check if the entry with name is managed
  inline bool contains( const std::string& name ) const { return theLayout.entries.count( name ); }

protected:
  /// An entry in the statevar vector consists of the name and a certain length
  struct StateVarEntryDefinition {
    std::string name;
    int         length;
  };

  /// The location in the statevar vector consists of the index and its certain length
  struct StateVarEntryLocation {
    int index;
    int length;
  };

  /// The layout is defined by a map of names to Locations, and the resulting required total length of the statevar
  /// vector
  struct StateVarVectorLayout {
    std::unordered_map< std::string, StateVarEntryLocation > entries;
    int                                                      nRequiredStateVars;
  };

  /// generate the statevar vector layout from a list of entries, defined by name and length
  static StateVarVectorLayout makeLayout( const std::vector< StateVarEntryDefinition >& theEntries )
  {
    std::unordered_map< std::string, StateVarEntryLocation > theMap;
    int                                                      sizeOccupied = 0;
    for ( const auto& theEntry : theEntries ) {
      const auto nextLocation = sizeOccupied;
      theMap[theEntry.name]   = { nextLocation, theEntry.length };
      sizeOccupied += theEntry.length;
    }
    return { theMap, sizeOccupied };
  }

  /// pointer to the first element in the statevar vector
  double* theStateVars;

  /// a const reference to the respective layout
  const StateVarVectorLayout& theLayout;

  MarmotStateVarVectorManager( double* theStateVars, const StateVarVectorLayout& theLayout_ )
    : theStateVars( theStateVars ), theLayout( theLayout_ ){};
};
