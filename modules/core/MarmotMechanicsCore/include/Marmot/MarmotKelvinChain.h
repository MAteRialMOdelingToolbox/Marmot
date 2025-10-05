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
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <functional>

namespace Marmot::Materials {

  namespace KelvinChain {
    /** This namespace contains the functions that derive the equivalent
     * Kelvin chain for a given compliance function
     * @typedef Properties
     * @brief Vector of material properties.
     *
     * Convenience typedef for an Eigen dynamic-size vector (`Eigen::VectorXd`)
     * that holds the material or model parameters of the Kelvin chain.
     */
    typedef Eigen::VectorXd Properties;

    /**
     * @typedef mapProperties
     * @brief Mapped view of material properties.
     *
     * Alias for `Eigen::Map<Properties>`, which allows mapping an existing
     * contiguous memory block as a `Properties` vector without copying.
     */
    typedef Eigen::Map< Properties > mapProperties;
    /**
     * @typedef StateVarMatrix
     * @brief Matrix of state variables.
     *
     * Alias for an Eigen matrix of shape 6 × N (`Eigen::Matrix<double, 6, Eigen::Dynamic>`),
     * storing viscoelastic strain variables for each Kelvin unit.
     */
    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;

    /**
     * @typedef mapStateVarMatrix
     * @brief Mapped view of the state variable matrix.
     *
     * Alias for `Eigen::Map<StateVarMatrix>`, allowing access to an existing
     * state variable array as an Eigen matrix without copying.
     */
    typedef Eigen::Map< StateVarMatrix > mapStateVarMatrix;
    /**
     * @brief Compile-time factorial.
     *
     * Recursive template structure computing the factorial of N at compile time.
     *
     * Example:
     * @code
     *   int f = Factorial<5>::value; // f = 120
     * @endcode
     *
     * @tparam N Non-negative integer whose factorial is to be computed.
     */
    template < int N >
    struct Factorial {
      enum { value = N * Factorial< N - 1 >::value };
    };

    template <>
    struct Factorial< 0 > {
      enum { value = 1 };
    };
    /**
     * Evaluates the Post-Widder formula for deriving the equivalent retardance \f$L_\kappa(\tau)\f$ of the discrete
     * Kelvin Chain model, given the continuous description of the compliance \f[\Phi(t) = \int_{\tau=0}^\infty
     * L(\tau)\left(1-e^{-\frac{t}{\tau}}\right)d(\ln\tau)\f] according to the formula
     * \f[L_k(\tau)=(-1)^{k-1}\frac{(-tk)^k\Phi^{(k)(t)}}{(k-1)!}\f]
     * in the sense that \f$L_k(\tau)=\lim_{\kappa\to\infty}L_\kappa(\tau)\f$ in the frequency domain.
     * @tparam k \f$k\f$ integer defining the order of differentiation that corresponds to the \f$L_k\f$ term in the
     * sequence of the retardance approximation
     * @param[in] phi compliance function whose discrete approximation is sought
     * @param[in] tau the continuous variable in the frequency domain
     * @return val the value of the \f$L_k(\tau)\f$ approximation.
     */
    template < int k >
    double evaluatePostWidderFormula( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                      double                                                                      tau )
    {
      autodiff::Real< k, double > tau_( tau * k );

      double val = -pow( -tau * k, k ) / double( Factorial< k - 1 >::value );
      val *= autodiff::derivatives( phi, autodiff::along( 1. ), autodiff::at( tau_ ) )[k];
      return val;
    }
    /**
     * Approximates the instantaneous elastic compliance \f$\Phi_0=\frac{1}{E_0}\f$,
     * \f[\frac{1}{E_0}=\int_{\tau=0}^{\tau}L(\tau)d\ln(\tau)\f]
     * It is used in order to avoid indegrating over the origin of the frequency domain where
     * \f$\int_0^{\tau_0}d\ln(\tau)\f$ is singular.
     * @tparam k \f$k\f$ integer defining the order of differentiation for the \f$k-th\f$ approximation of the
     * retardance.
     * @param[in] phi compliance function whose discrete approximation is sought.
     * @param[in] tauMin the minimum retardation time that is required to be of the same order or smaller than the total
     * time of the analysis.
     * @param[in] spacing the increment between two sunsequent retardation times
     * @return val the value of the zeroth compliance term.
     */

    template < int k >
    double approximateZerothCompliance( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                        double tauMin,
                                        double spacing = 10. )
    {
      NumericalAlgorithms::Integration::scalar_to_scalar_function_type f = [&]( double tau ) {
        double                      val_ = -pow( -k, k ) * pow( tau, k - 1 ) / double( Factorial< k - 1 >::value );
        autodiff::Real< k, double > tau_( tau * k );
        val_ *= autodiff::derivatives( phi, autodiff::along( 1. ), autodiff::at( tau_ ) )[k];
        return val_;
      };

      double
        val = NumericalAlgorithms::Integration::integrateScalarFunction( f,
                                                                         { 1e-14, tauMin / sqrt( spacing ) },
                                                                         100,
                                                                         NumericalAlgorithms::Integration::simpson );
      return val;
    }
    /**
     * Evaluates the elastic moduli based on the \f$L_k\f$ approximation of the retardance,
     * \f[ \frac{1}{E_\mu}=(\ln 10)L_k(\tau_\mu),\;\mu=1,2,...,M\f]
     * It is used in order to avoid indegrating over the origin of the frequency domain where
     * \f$\int_0^{\tau_0}d\ln(\tau)\f$ is singular.
     * @tparam k \f$k\f$ integer defining the order of differentiation for the \f$k-th\f$ approximation of the
     * retardance.
     * @param[in] phi compliance function whose discrete approximation is sought
     * @param[in] retardationTimes the minimum retardation time that is required to be of the same order or smaller than
     * the total time of the analysis
     * @param[in] gaussQuadrature flag if on Gauss Quadrature is performed else the Post-Widder formula is applied.
     * @returns elasticModuli the discrete elastic moduli of the equivalt Kelvin chain.
     */

    template < int k >
    Properties computeElasticModuli( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                     Properties retardationTimes,
                                     bool       gaussQuadrature = false )
    {
      Properties elasticModuli( retardationTimes.size() );
      double     spacing = retardationTimes( 1 ) / retardationTimes( 0 );

      for ( int i = 0; i < retardationTimes.size(); i++ ) {
        double tau = retardationTimes( i );
        if ( !gaussQuadrature ) {
          elasticModuli( i ) = 1. / ( log( spacing ) * evaluatePostWidderFormula< k >( phi, tau ) );
        }
        else {
          elasticModuli( i ) = 1. /
                               ( log( spacing ) / 2. *
                                 ( evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. ) / 6. ) ) +
                                   evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. ) / 6. ) ) ) );
        }
      }

      return elasticModuli;
    }
    /**
     * This function evaluates the retardation times.
     * @param[in] n the number of kelvin units in the equivalent Kelvin chain.
     * @param[in] min the first retardation time.
     * @param[in] spacing the spacing between the retardation times.
     * @returns retardationTimes a vector containing the retardation time for each Kelvin unit in the Kelvin chain.
     */

    Properties generateRetardationTimes( int n, double min, double spacing );
    /**
     * For the given time increment \f$\Delta t_k\f$ this function updates the visco-elastic state variables according
     * to the update rule: \f[\varepsilon^{ev,\mu,k+1}_{ij} = \frac{\lambda^{\mu
     * k}}{E_\mu}C^\nu_{ijkl}\Delta\sigma_{kl}+\beta^{\mu ,k}\varepsilon^{ev,\mu, k}_{ij}\f]
     * @param[in] dT the time increment.
     * @param[in] elasticModuli vector containing the elastic moduli of ech Kelvin unit in the Kelvin chain.
     * @param[in] retardationTimes vector containing the retardation time for each Kelvin unit in the Kelvin chain.
     * @param[in,out] stateVars the \f$[6\times \mu]\f$ matrix that contains the viscoelastic strain update for each
     * unit of the Kelvin chain.
     * @param[in] dStress the \f$[6\times 1]\f$ vector of the total stress increment.
     * @param[in] unitComplianceMatrix the [6\times 6] compliance matrix of the material with unit compliance and given
     * Poisson coeffiscient.
     */

    void updateStateVarMatrix( const double                 dT,
                               Properties                   elasticModuli,
                               Properties                   retardationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStress,
                               const Marmot::Matrix6d&      unitComplianceMatrix );
    /**
     * For the given time increment \f$\Delta t_k\f$ this function updates the visco-elastic state variables according
     *to the update rule: \f[\overline{J}^k=\sum_{\mu=1}^M\frac{1-\lambda^{\mu,k}}{E^\mu}\f]
     *\f[\Delta\varepsilon^{'',k}_{ij} = \sum_{\mu=1}^M\left(1-\beta^{\mu,k}\right)\varepsilon^{ev,\kappa}_{ij}\f]
     * @param[in] dT the time increment.
     * @param[in] elasticModuli vector containing the elastic moduli of ech Kelvin unit in the Kelvin chain.
     * @param[in] retardationTimes vector containing the retardation time for each Kelvin unit in the Kelvin chain.
     * @param[in] stateVars the \f$[6\times \mu]\f$ matrix that contains the viscoelastic strain update for each unit of
     *the Kelvin chain.
     * @param[in,out] uniaxialCompliance number containing the \f$\sum^M_{\mu=1}\frac{\lambda^{\mu, k}}{E_\mu}\f$.
     * @param[in,out] dStrain the viscoelastic update of the strain based on the state variables of the Kelvin chain.
     * @param[in] factor the solidification factor (in non aging viscoelasticity set to 1).
     */

    void evaluateKelvinChain( const double      dT,
                              Properties        elasticModuli,
                              Properties        retardationTimes,
                              StateVarMatrix    stateVars,
                              double&           uniaxialCompliance,
                              Marmot::Vector6d& dStrain,
                              const double      factor );
    /**
     * This function evaluates the time parameters for each Kelvin chain.
     * @param[in] dT the time increment.
     * @param[in] tau the retardation time for a given Kelvin unit.
     * @param[out] lambda time dependent factor for each Kelvin unit.
     * @param[out] beta time dependemnt factor for each Kelvin unit.
     */

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace KelvinChain
} // namespace Marmot::Materials
