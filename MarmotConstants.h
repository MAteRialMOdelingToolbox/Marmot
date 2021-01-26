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

namespace Marmot {
    namespace Constants {
        constexpr double    Pi         = 3.141592653589793238463;
        constexpr double    numZeroPos = 1e-16;
        double              cubicRootEps();
        double              squareRootEps();
        extern const double CubicRootEps;
        extern const double SquareRootEps;
        constexpr double    sqrt3_8 = 0.61237243569579452454932101867647;
        constexpr double    sqrt2_3 = 0.8164965809277260327324280249019;
        constexpr double    sqrt3_2 = 1.2247448713915890490986420373529;
        constexpr double    sqrt2   = 1.4142135623730950488016887242097;
        constexpr double    sqrt3   = 1.7320508075688772935274463415059;
        constexpr double    sqrt6   = 2.4494897427831780981972840747059;
    } // namespace Constants
} // namespace Marmot
