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
#include <cmath>
#include <tuple>

namespace Marmot::GradientDamage {

  namespace DecreasingInteractions {

    /**
     * @brief Damage dependent nonlocal interaction function acc. Poh & Sun (2017)
     * @param omega damage variable
     * @param eta parameter controlling the rate of decrease of interactions with damage
     * @param R parameter controlling the residual interactions at full damage
     * @tparam T type of the interaction function value (e.g., double, dual, etc.)
     * @return interaction function value
     *
     * @details The function is defined as:
     * \f[ g(\omega) = \frac{(1 - R) \exp(-\eta \, \omega) + R - \exp(-\eta) }{ 1 - \exp(-\eta)} \f]
     * where \f$ \omega \f$ is the damage variable, \f$ \eta \f$ controls the rate of decrease of interactions with
     * damage, and \f$ R \f$ controls the residual interactions at full damage. The function is designed to decrease as
     * damage increases, with a residual interaction value of \f$ R \f$ at full damage (\f$ \omega = 1 \f$).
     *
     */
    template < typename T >
    T poh( T omega, double eta, double R )
    {
      const T g = ( ( 1. - R ) * exp( -eta * omega ) + R - exp( -eta ) ) / ( 1. - exp( -eta ) );
      return g;
    }

    namespace FirstOrderDerived {

      /**
       * @brief First-order derivative of the Poh & Sun (2017) interaction function with respect to damage
       * @param omega damage variable
       * @param eta parameter controlling the rate of decrease of interactions with damage
       * @param R parameter controlling the residual interactions at full damage
       * @return a tuple containing the interaction function value and its first derivative with respect to damage
       * @details The first-order derivative variant of the interaction function @p poh.
       */
      std::tuple< double, double > poh( double omega, double eta, double R );

    } // namespace FirstOrderDerived

    namespace SecondOrderDerived {

      /**
       * @brief Second-order derivative of the Poh & Sun (2017) interaction function with respect to damage
       * @param omega damage variable
       * @param eta parameter controlling the rate of decrease of interactions with damage
       * @param R parameter controlling the residual interactions at full damage
       * @return a tuple containing the interaction function value, its first derivative, and its second derivative with
       * respect to damage
       * @details The second-order derivative variant of the interaction function @p poh.
       */
      std::tuple< double, double, double > poh( double omega, double eta, double R );

    } // namespace SecondOrderDerived
  }   // namespace DecreasingInteractions

} // namespace Marmot::GradientDamage
