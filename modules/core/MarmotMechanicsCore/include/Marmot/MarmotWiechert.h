/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ ___ ___   ___ | |_
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

#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"

#include <cmath>
#include <functional>

namespace Marmot::Materials {

  namespace Wiechert {

    using Properties        = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::Properties;
    using mapProperties     = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::mapProperties;
    using StateVarMatrix    = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::StateVarMatrix;
    using mapStateVarMatrix = Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::mapStateVarMatrix;

    /**
     * @brief Evaluate the Post-Widder formula of order @f$k@f$ for a relaxation function.
     *
     * Approximates the relaxation spectrum of the relaxation function @f$\psi@f$ at the relaxation time @f$\tau@f$.
     * The derivative of order @f$k@f$ is obtained by automatic differentiation.
     *
     * @tparam k order of the Post-Widder approximation
     * @param[in] psi relaxation function
     * @param[in] tau relaxation time at which the spectrum is evaluated
     * @return the relaxation spectrum at @f$\tau@f$
     */
    template < int k >
    double evaluatePostWidderFormula( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > psi,
                                      double                                                                      tau )
    {
      using Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::PostWidderCoefficientSign;
      return Marmot::ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::evaluatePostWidderFormula<
        k >( psi, tau, PostWidderCoefficientSign::Positive );
    }

    /**
     * @brief Compute the elastic moduli of the Maxwell units of a Wiechert chain.
     *
     * The relaxation spectrum is integrated over the logarithmic relaxation-time axis, which is discretized by the
     * relaxation times with the given spacing. By default a one-point rule is used, i.e. the modulus of the unit with
     * relaxation time @f$\tau_i@f$ is @f$\ln(\mathrm{spacing}) \, S(\tau_i)@f$, where @f$S@f$ is the relaxation
     * spectrum. With Gauss quadrature the two-point Gauss rule over the logarithmic interval of the unit is used.
     *
     * @tparam k order of the Post-Widder approximation of the relaxation spectrum
     * @param[in] psi relaxation function
     * @param[in] relaxationTimes relaxation times of the Maxwell units
     * @param[in] spacing ratio of two consecutive relaxation times
     * @param[in] gaussQuadrature use the two-point Gauss rule instead of the one-point rule
     * @return the elastic moduli, one per Maxwell unit
     */
    template < int k >
    Properties computeElasticModuli( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > psi,
                                     const Properties& relaxationTimes,
                                     double            spacing,
                                     bool              gaussQuadrature = false )
    {
      Properties elasticModuli( relaxationTimes.size() );
      for ( int i = 0; i < relaxationTimes.size(); ++i ) {
        const double tau = relaxationTimes( i );
        if ( !gaussQuadrature ) {
          elasticModuli( i ) = std::log( spacing ) * evaluatePostWidderFormula< k >( psi, tau );
        }
        else {
          elasticModuli(
            i ) = std::log( spacing ) / 2. *
                  ( evaluatePostWidderFormula< k >( psi, tau * std::pow( spacing, -std::sqrt( 3. ) / 6. ) ) +
                    evaluatePostWidderFormula< k >( psi, tau * std::pow( spacing, std::sqrt( 3. ) / 6. ) ) );
        }
      }
      return elasticModuli;
    }

    /**
     * @brief Generate logarithmically spaced relaxation times.
     *
     * The relaxation times are @f$\tau_i = \mathrm{min} \cdot \mathrm{spacing}^i@f$ for @f$i = 0, \dots, n-1@f$.
     *
     * @param[in] n number of relaxation times
     * @param[in] min smallest relaxation time
     * @param[in] spacing ratio of two consecutive relaxation times
     * @return the relaxation times
     */
    Properties generateRelaxationTimes( int n, double min, double spacing );

    /**
     * @brief Update the stresses of the Maxwell units for a strain increment.
     *
     * For every unit @f$i@f$ the stress is updated with the exponential algorithm,
     * @f$\sigma_i \leftarrow \lambda_i D_i \, \mathbb{C}_1 \, \Delta\varepsilon + \beta_i \sigma_i@f$,
     * where @f$\lambda_i@f$ and @f$\beta_i@f$ follow from the time increment and the relaxation time of the unit.
     * Nothing is updated for a vanishing time increment.
     *
     * @param[in] dT time increment
     * @param[in] elasticModuli elastic moduli @f$D_i@f$ of the Maxwell units
     * @param[in] relaxationTimes relaxation times of the Maxwell units
     * @param[in,out] stateVars stresses of the Maxwell units in Voigt notation, one column per unit
     * @param[in] dStrain strain increment in Voigt notation
     * @param[in] unitD_ijkl stiffness tensor @f$\mathbb{C}_1@f$ for a unit Young's modulus in Voigt notation
     */
    void updateStateVarMatrix( const double                 dT,
                               const Properties&            elasticModuli,
                               const Properties&            relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStrain,
                               const Marmot::Matrix6d&      unitD_ijkl );

    /**
     * @brief Accumulate the contribution of the Maxwell units to the stiffness and to the stress relaxation.
     *
     * For every unit @f$i@f$ the algorithmic stiffness @f$\lambda_i D_i@f$ is added to @p uniaxialStiffness and the
     * stress @f$(1 - \beta_i) \sigma_i@f$ released by relaxation during the increment is added to @p dStress,
     * both multiplied by @p factor.
     *
     * @param[in] dT time increment
     * @param[in] elasticModuli elastic moduli @f$D_i@f$ of the Maxwell units
     * @param[in] relaxationTimes relaxation times of the Maxwell units
     * @param[in] stateVars stresses @f$\sigma_i@f$ of the Maxwell units at the beginning of the increment
     * @param[in,out] uniaxialStiffness accumulated algorithmic stiffness of the Maxwell units
     * @param[in,out] dStress accumulated stress released by relaxation, in Voigt notation
     * @param[in] factor scalar weight of the contribution
     */
    void evaluateWiechert( const double                 dT,
                           const Properties&            elasticModuli,
                           const Properties&            relaxationTimes,
                           Eigen::Ref< StateVarMatrix > stateVars,
                           double&                      uniaxialStiffness,
                           Marmot::Vector6d&            dStress,
                           const double                 factor );

  } // namespace Wiechert
} // namespace Marmot::Materials
