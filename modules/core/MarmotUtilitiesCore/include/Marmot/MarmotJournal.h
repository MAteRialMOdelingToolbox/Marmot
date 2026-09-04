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
#include <iostream>
#include <sstream>
#include <string>

namespace Marmot {

  /** @class MakeString
   * @brief Utility class for constructing strings with stream-like syntax.
   *
   * This class allows for easy string construction using the stream insertion operator.
   * It can be used to build complex strings in a readable manner.
   */
  class MakeString {
  public:
    std::stringstream stream; ///< Internal string stream used for building the output string.

    /// Convert the accumulated stream content to a `std::string`.
    operator std::string() const { return stream.str(); }

    /**
     * @brief Append a value to the internal stream.
     * @tparam T   Type of the value to append.
     * @param VAR  Value to insert into the stream.
     * @return     Reference to `*this` for chaining.
     */
    template < class T >
    MakeString& operator<<( T const& VAR )
    {
      stream << VAR;
      return *this;
    }
  };

  /**
   * @class MarmotJournal
   * @brief Singleton class for managing output messages in the Marmot framework.
   *
   * This class provides a centralized way to handle warning and notification messages,
   * allowing them to be directed to a specified output stream (e.g., console, file).
   */
  class MarmotJournal {
  private:
    static MarmotJournal& getInstance();

    std::ostream output;

    MarmotJournal();

  public:
    MarmotJournal( MarmotJournal const& )  = delete;
    void operator=( MarmotJournal const& ) = delete;

    /**
     * @brief Redirect all subsequent journal output to @p newOutputStream.
     * @param newOutputStream  The output stream that warnings and notifications are written to.
     */
    static void setMSGOutputDirection( std::ostream& newOutputStream );

    /**
     * @brief Write a warning message to the journal output stream.
     * @param message  The warning text to emit.
     * @return `true` after the message has been written.
     */
    static bool warningToMSG( const std::string& message );

    /**
     * @brief Write a notification message to the journal output stream.
     * @param message  The notification text to emit.
     * @return `true` after the message has been written.
     */
    static bool notificationToMSG( const std::string& message );
  };

} // namespace Marmot
