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
#include <cmath>

namespace Marmot {
  namespace Constants {
    constexpr double Pi         = 3.141592653589793238463;
    constexpr double numZeroPos = 1e-16;

    /** Inlined function to get the cubic root of machine epsilon for double precision.
     *
     * This function computes and returns the cubic root of the machine epsilon for double
     * precision floating-point numbers. Machine epsilon is the smallest positive number that,
     * when added to 1.0, results in a value different from 1.0 due to rounding errors.
     *
     * @return The cubic root of machine epsilon for double precision.
     */
    inline double cubicRootEps()
    {
      return std::pow( std::numeric_limits< double >::epsilon(), 1. / 3 );
    }

    /** @brief Inlined function to get the square root of machine epsilon for double precision.
     *
     *  This function computes and returns the square root of the machine epsilon for double
     *  precision floating-point numbers. Machine epsilon is the smallest positive number that,
     *  when added to 1.0, results in a value different from 1.0 due to rounding errors.
     *
     *  @return The square root of machine epsilon for double precision.
     */
    inline double squareRootEps()
    {
      return std::pow( std::numeric_limits< double >::epsilon(), 0.5 );
    }

    /** The golden ratio, approximately 1.618033988749895.
     *
     * The golden ratio is a mathematical constant that appears in various areas of mathematics,
     * art, and nature. It is defined as (1 + sqrt(5)) / 2 and is often denoted by the Greek letter phi (φ).
     */
    constexpr double GoldenRatio = 1.618033988749895;

    /** The golden angle in radians, approximately 2.399963229728653.
     *
     * The golden angle is the angle that divides a circle in such a way that the ratio of the
     * larger arc to the smaller arc is equal to the golden ratio. It is approximately 137.5 degrees
     * or 2.399963229728653 radians.
     */
    constexpr double GoldenAngle = 2.399963229728653;

    const inline double SquareRootEps = squareRootEps();
    const inline double CubicRootEps  = cubicRootEps();

    /** \@brief Constant expression for \f$\sqrt{3/8}$ */
    constexpr double sqrt3_8 = 0.61237243569579452454932101867647;

    /** \@brief Constant expression for \f$\sqrt{2/3}$ */
    constexpr double sqrt2_3 = 0.8164965809277260327324280249019;
    constexpr double sqrt3_2 = 1.2247448713915890490986420373529;
    constexpr double sqrt2   = 1.4142135623730950488016887242097;
    constexpr double sqrt3   = 1.7320508075688772935274463415059;
    constexpr double sqrt6   = 2.4494897427831780981972840747059;
  } // namespace Constants
} // namespace Marmot
