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

namespace Marmot::PhaseField {

  namespace EnergyDegradationFunctions {
    // degradation functions acc. Kuhn, Schlüter, Müller (2015)

    /**
     * @brief Quadratic degradation function.
     * @param pf Phase field variable.
     * @tparam T Numeric type (e.g., double, float).
     * @return Degradation function value.
     * @details This function is defined as
     * \f[ g(\phi) = (1 - \phi)^2,\f]
     * where \f$\phi\f$ is the phase field variable.
     * It represents a quadratic degradation of the material stiffness as the phase field variable increases.
     */
    template < typename T >
    T quadratic( const T pf )
    {
      return pow( 1. - pf, 2 );
    }

    /**
     * @brief Cubic degradation function.
     * @param pf Phase field variable.
     * @tparam T Numeric type (e.g., double, float).
     * @return Degradation function value.
     * @details This function is defined as
     * \f[ g(\phi) = 3(1 - \phi)^2 - 2(1 - \phi)^3,\f]
     * where \f$\phi\f$ is the phase field variable.
     * It represents a cubic degradation of the material stiffness as the phase field variable increases,
     * providing a smoother transition compared to the quadratic function.
     */
    template < typename T >
    T cubic( const T pf )
    {
      const T s = 1. - pf;
      return 3. * pow( s, 2 ) - 2. * pow( s, 3 );
    }

    /**
     * @brief Quartic degradation function.
     * @param pf Phase field variable.
     * @tparam T Numeric type (e.g., double, float).
     * @return Degradation function value.
     * @details This function is defined as
     * \f[ g(\phi) = 4(1 - \phi)^3 - 3(1 - \phi)^4,\f]
     * where \f$\phi\f$ is the phase field variable.
     * It represents a quartic degradation of the material stiffness as the phase field variable increases,
     * providing an even smoother transition compared to the cubic function.
     */
    template < typename T >
    T quartic( const T pf )
    {
      const T s = 1. - pf;
      return 4 * pow( s, 3 ) - 3 * pow( s, 4 );
    }

    /**
     * @brief Generic degradation function acc. to Wu, JMPS (2017).
     * @param pf Phase field variable.
     * @param p Exponent parameter controlling the shape of the degradation function.
     * @param a1 Coefficient for the linear term in the denominator.
     * @param a2 Coefficient for the quadratic term in the denominator.
     * @param a3 Coefficient for the cubic term in the denominator.
     * @tparam T Numeric type (e.g., double, float).
     * @return Degradation function value.
     * @details This function is defined as
     * \f[ g(\phi) = \frac{(1 - \phi)^p}{(1 - \phi)^p + a_1 \phi (1 + a_2 \phi (1 + a_3 \phi))},\f]
     * where \f$\phi\f$ is the phase field variable, and \f$p\f$, \f$a_1\f$, \f$a_2\f$, and \f$a_3\f$ are parameters
     * that control the shape of the degradation function. This generic form allows for a wide range of degradation
     * behaviors by adjusting the parameters accordingly.
     */
    template < typename T >
    T generic( const T pf, const double p, const double a1, const double a2, const double a3 )
    {
      const T s      = 1. - pf;
      const T P      = 1. + a2 * pf * ( 1. + a3 * pf );
      const T Q      = a1 * pf * P;
      const T result = pow( s, p ) / ( pow( s, p ) + Q );

      return result;
    }

    namespace SecondOrderDerived {

      /**
       * @brief Quadratic degradation function and its first and second derivatives.
       * @param pf Phase field variable.
       * @return A tuple containing the degradation function value, its first derivative, and its second derivative with
       * respect to the phase field variable.
       * @details This function returns the value of the quadratic degradation function along with its first and second
       * derivatives. The degradation function is defined as \f[ g(\phi) = (1 - \phi)^2,\f] and its derivatives are
       * given by \f[ g'(\phi) = -2(1 - \phi), \quad g''(\phi) = 2. \f]
       */
      std::tuple< double, double, double > quadratic( const double pf );

      /**
       * @brief Cubic degradation function and its first and second derivatives.
       * @param pf Phase field variable.
       * @return A tuple containing the degradation function value, its first derivative, and its second derivative with
       * respect to the phase field variable.
       * @details This function returns the value of the cubic degradation function along with its first and second
       * derivatives. The degradation function is defined as \f[ g(\phi) = 3(1 - \phi)^2 - 2(1 - \phi)^3,\f] and its
       * derivatives are given by \f[ g'(\phi) = -6(1 - \phi)\phi, \quad g''(\phi) = 12\phi - 6. \f]
       */
      std::tuple< double, double, double > cubic( const double pf );

      /**
       * @brief Quartic degradation function and its first and second derivatives.
       * @param pf Phase field variable.
       * @return A tuple containing the degradation function value, its first derivative, and its second derivative with
       * respect to the phase field variable.
       * @details This function returns the value of the quartic degradation function along with its first and second
       * derivatives. The degradation function is defined as \f[ g(\phi) = 4(1 - \phi)^3 - 3(1 - \phi)^4,\f] and its
       * derivatives are given by \f[ g'(\phi) = -12(1 - \phi)^2 + 12(1 - \phi)^3, \quad g''(\phi) = 24(1 - \phi) -
       * 36(1 - \phi)^2. \f]
       */
      std::tuple< double, double, double > quartic( const double pf );

      /**
       * @brief Generic degradation function and its first and second derivatives.
       * @param pf Phase field variable.
       * @param p Exponent parameter controlling the shape of the degradation function.
       * @param a1 Coefficient for the linear term in the denominator.
       * @param a2 Coefficient for the quadratic term in the denominator.
       * @param a3 Coefficient for the cubic term in the denominator.
       * @return A tuple containing the degradation function value, its first derivative, and its second derivative with
       * respect to the phase field variable.
       * @details This function returns the value of the generic degradation function along with its first and second
       * derivatives. The degradation function is defined as \f[ g(\phi) = \frac{(1 - \phi)^p}{(1 - \phi)^p + a_1 \phi
       * (1 + a_2 \phi (1 + a_3 \phi))},\f] where \f$\phi\f$ is the phase field variable, and \f$p\f$, \f$a_1\f$,
       * \f$a_2\f$, and \f$a_3\f$ are parameters that control the shape of the degradation function. The first and
       * second derivatives are computed using automatic differentiation with hyper-dual numbers.
       */
      std::tuple< double, double, double > generic( const double pf,
                                                    const double p,
                                                    const double a1,
                                                    const double a2,
                                                    const double a3 );
    }; // namespace SecondOrderDerived
  };   // namespace EnergyDegradationFunctions
};     // namespace Marmot::PhaseField
