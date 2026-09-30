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

/** @file MarmotPortability.h
 * @brief Compiler portability shims, included by the Marmot utility headers that (almost) every Marmot file uses.
 *
 * @details MSVC does not provide the GCC/Clang extension `__PRETTY_FUNCTION__`, which Marmot uses throughout its
 * error messages. MSVC's equivalent is `__FUNCSIG__`.
 */

#if defined( _MSC_VER ) && !defined( __PRETTY_FUNCTION__ )
/** @brief Full signature of the enclosing function, mapped to MSVC's equivalent of the GCC/Clang extension. */
#  define __PRETTY_FUNCTION__ __FUNCSIG__
#endif
