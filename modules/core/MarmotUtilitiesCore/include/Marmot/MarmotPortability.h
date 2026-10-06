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
 * error messages. MSVC's equivalent is `__FUNCSIG__`. `__PRETTY_FUNCTION__` is a predefined identifier, not a macro,
 * so its presence cannot be tested with `defined`; the shim is keyed on the compiler instead, and excludes clang-cl,
 * which defines `_MSC_VER` but provides `__PRETTY_FUNCTION__` natively.
 *
 * MARMOT_API marks what the Marmot library exports to its consumers (EdelweissFE, the Abaqus/CADFEM interfaces).
 * A Windows DLL exports only what is marked, and only that is available to code linking Marmot; everything else
 * is reached through virtual functions of objects created by the Marmot factories. On other platforms, everything
 * is exported unless Marmot is built with MARMOT_EXPORT_API_ONLY, which mimics the Windows behavior.
 */

#if defined( _WIN32 )
#  if defined( MARMOT_BUILDING_LIBRARY )
/** @brief Exports a class or function from the Marmot library (defined while building it). */
#    define MARMOT_API __declspec( dllexport )
#  else
/** @brief Imports a class or function from the Marmot library. */
#    define MARMOT_API __declspec( dllimport )
#  endif
#else
/** @brief Exports a class or function from the Marmot library, also when it is built with hidden visibility. */
#  define MARMOT_API __attribute__( ( visibility( "default" ) ) )
#endif

#if defined( _MSC_VER ) && !defined( __clang__ )
/** @brief Full signature of the enclosing function, mapped to MSVC's equivalent of the GCC/Clang extension. */
#  define __PRETTY_FUNCTION__ __FUNCSIG__
#endif
