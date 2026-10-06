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
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "autodiff/forward/real.hpp"

#include <cmath>
#include <functional>

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

      namespace RelaxationFunctions {

        /**
         * @brief Power-law relaxation function.
         *
         * Computes the relaxation function
         * \f[
         *   \Psi(\tau) = m \, \tau^{-n}.
         * \f]
         *
         * @tparam T_ Scalar or autodiff type of the argument.
         * @param[in] tau Relaxation time or evaluation point.
         * @param[in] m Scaling factor.
         * @param[in] n Positive exponent controlling the relaxation rate.
         * @return The relaxation-function value \f$\Psi(\tau)\f$.
         */
        template < typename T_ >
        T_ powerLaw( T_ tau, double m, double n )
        {
          return m * pow( tau, -n );
        }

      } // namespace RelaxationFunctions

      namespace DiscreteSpectrum {

        using Properties        = Eigen::VectorXd;
        using mapProperties     = Eigen::Map< Properties >;
        using StateVarMatrix    = Eigen::Matrix< double, 6, Eigen::Dynamic >;
        using mapStateVarMatrix = Eigen::Map< StateVarMatrix >;

        /// Sign of the coefficient of the Post-Widder formula: positive for relaxation, negative for retardation.
        enum class PostWidderCoefficientSign { Positive, Negative };

        /**
         * @brief Evaluate the Post-Widder formula of order @f$k@f$.
         *
         * Computes @f$\pm \frac{(-k\tau)^k}{(k-1)!} \, f^{(k)}(k\tau)@f$, with the derivative of order @f$k@f$
         * obtained by automatic differentiation.
         *
         * @tparam k order of the approximation
         * @param[in] function function to be inverted, e.g. a relaxation or retardation function
         * @param[in] tau relaxation or retardation time at which the spectrum is evaluated
         * @param[in] coefficientSign sign of the coefficient
         * @return the value of the approximated spectrum at @f$\tau@f$
         */
        template < int k >
        double evaluatePostWidderFormula(
          std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > function,
          double                                                                      tau,
          PostWidderCoefficientSign                                                   coefficientSign )
        {
          autodiff::Real< k, double > evaluationTime( tau * k );
          double coefficient = std::pow( -tau * k, k ) / static_cast< double >( Marmot::Math::factorial( k - 1 ) );
          if ( coefficientSign == PostWidderCoefficientSign::Negative )
            coefficient *= -1.;
          return coefficient *
                 autodiff::derivatives( function, autodiff::along( 1. ), autodiff::at( evaluationTime ) )[k];
        }

        /**
         * @brief Generate logarithmically spaced relaxation or retardation times.
         *
         * The times are @f$\tau_i = \mathrm{min} \cdot \mathrm{spacing}^i@f$ for @f$i = 0, \dots, n-1@f$.
         *
         * @param[in] n number of times
         * @param[in] min smallest time
         * @param[in] spacing ratio of two consecutive times
         * @return the times
         */
        Properties generateLogarithmicTimes( int n, double min, double spacing );

        /**
         * @brief Compute the integration factors of the exponential algorithm for a Maxwell or Kelvin unit.
         *
         * With @f$x = \Delta t / \tau@f$ the factors are @f$\beta = e^{-x}@f$ and
         * @f$\lambda = (1 - \beta) / x@f$. For very small and very large @f$x@f$ the limiting expressions are used
         * (Jirasek and Bazant).
         *
         * @param[in] dT time increment
         * @param[in] tau relaxation or retardation time of the unit
         * @param[out] lambda factor @f$\lambda@f$
         * @param[out] beta factor @f$\beta@f$
         */
        void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

      } // namespace DiscreteSpectrum
    }   // namespace Viscoelasticity
  }     // namespace ContinuumMechanics
} // namespace Marmot
