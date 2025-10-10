
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
#include "Marmot/MarmotVoigt.h"

namespace Marmot {

  namespace ContinuumMechanics {

    namespace Viscoelasticity {

      namespace ComplianceFunctions {

        /**
         * @brief Logarithmic power-law compliance function.
         *
         * Computes the compliance function
         * \f[
         *   \Phi(\tau) = m \, \ln\!\left( 1 + \tau^n \right),
         * \f]
         * creep for long times.
         *
         * @tparam T_ Scalar or autodiff type of the argument.
         *
         * @param[in] tau Retardation time or evaluation point.
         * @param[in] m Scaling factor.
         * @param[in] n Exponent controlling the growth rate.
         *
         * @return The compliance value \f$\Phi(\tau)\f$.
         */

        template < typename T_ >
        T_ logPowerLaw( T_ tau, double m, double n )
        {
          T_ val = m * log( 1. + pow( tau, n ) );
          return val;
        }
        /**
         * @brief Power-law compliance function.
         *
         * Computes the compliance function
         * \f[
         *   \Phi(\tau) = m \, \tau^n,
         * \f]
         * creep for short times.
         *
         * @tparam T_ Scalar or autodiff type of the argument.
         *
         * @param[in] tau Retardation time or evaluation point.
         * @param[in] m Scaling factor.
         * @param[in] n Exponent controlling the growth rate.
         *
         * @return The compliance value \f$\Phi(\tau)\f$.
         */
        template < typename T_ >
        T_ powerLaw( T_ tau, double m, double n )
        {
          T_ val = m * pow( tau, n );
          return val;
        }
      } // namespace ComplianceFunctions
    }   // namespace Viscoelasticity
  }     // namespace ContinuumMechanics
} // namespace Marmot
