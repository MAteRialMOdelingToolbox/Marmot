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

/*
 * Fastor/Fastor.h ends with `#pragma warning( default : 4100 )` (among others), which re-enables the
 * unused-parameter warning for everything that follows, and so overrides the /wd4100 given on the command line.
 * Unused parameters are intended in the interface implementations (cf. -Wno-unused-parameter for GCC/Clang), so
 * this header has to be included right after Fastor/Fastor.h.
 */
#ifdef _MSC_VER
#  pragma warning( disable : 4100 )
#endif
