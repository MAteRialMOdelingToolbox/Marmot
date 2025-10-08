
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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include <functional>

namespace Marmot {
  namespace NumericalAlgorithms::Integration {

    /** @brief Type alias for a function that takes a double and returns a double.*/
    using scalar_to_scalar_function_type = std::function< double( const double x ) >;

    /** @brief Enumeration of available numerical integration rules.
     *
     *  This enum defines the different numerical integration methods that can be used
     *  in the integrateScalarFunction function.
     */
    enum integrationRule {
      /// Midpoint rule for numerical integration
      midpoint,
      /// Trapezoidal rule for numerical integration
      trapezodial,
      /// Simpson's rule for numerical integration
      simpson
    };

    /** @brief Numerically integrates a scalar function over a specified interval using a chosen integration rule.
     *  @param f The scalar function to be integrated.
     *  @param integrationLimits A tuple containing the lower and upper limits of integration.
     *  @param n The number of subintervals to use for the integration (must be positive).
     *  @param intRule The numerical integration rule to use (midpoint, trapezoidal, or Simpson's rule).
     *  @return The approximate value of the integral of f over the specified interval.
     *
     *  This function approximates the integral of the given scalar function @p f over the interval
     *  defined by integrationLimits using the specified numerical integration rule. The number
     *  of subintervals n determines the accuracy of the approximation; a larger @p n generally
     *  leads to a more accurate result.
     *
     *  Example usage:
     *  ```cpp
     *  auto f = [](double x) { return x * x; }; // Function to integrate
     *  double result = integrateScalarFunction(f, std::make_tuple(0.0, 1.0), 100, integrationRule::simpson);
     *  ```
     */
    double integrateScalarFunction( scalar_to_scalar_function_type     f,
                                    const std::tuple< double, double > integrationLimits,
                                    const int                          n,
                                    const integrationRule              intRule );
  } // namespace NumericalAlgorithms::Integration
} // namespace Marmot
