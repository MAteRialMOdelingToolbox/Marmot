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
#include <iostream>
#include <sstream>
#include <string>

class MakeString {
public:
  std::stringstream stream;
                    operator std::string() const { return stream.str(); }

  template < class T >
  MakeString& operator<<( T const& VAR )
  {
    stream << VAR;
    return *this;
  }
};

class MarmotJournal {
private:
  static MarmotJournal& getInstance();

  std::ostream output;

  MarmotJournal();

public:
  MarmotJournal( MarmotJournal const& )  = delete;
  void operator=( MarmotJournal const& ) = delete;

  static void setMSGOutputDirection( std::ostream& newOutputStream );

  static bool warningToMSG( const std::string& message );

  static bool notificationToMSG( const std::string& message );
};
