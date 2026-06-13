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
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
#include "autodiff/forward/real.hpp"

#include <functional>

namespace Marmot::Materials {

  /**
   * @namespace Wiechert
   * @brief Utilities for branch-wise Wiechert viscoelastic updates in 3D Voigt notation.
   *
   * @details
   * This namespace provides helper functions to initialize Maxwell branch properties,
   * update branch state variables, evaluate incremental viscoelastic response contributions,
   * and compute robust integration factors \f$\lambda\f$ and \f$\beta\f$.
   */
  namespace Wiechert {

    /// @brief Dynamic vector of branch-wise material properties (e.g. moduli or relaxation times).
    typedef Eigen::VectorXd Properties;

    /// @brief Non-owning Eigen map view onto a `Properties` buffer.
    typedef Eigen::Map< Properties > mapProperties;

    /// @brief Branch state variable matrix (6 x nMaxwell) storing Voigt stress-like internal variables.
    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;

    /// @brief Non-owning Eigen map view of `StateVarMatrix`.
    typedef Eigen::Map< StateVarMatrix > mapStateVarMatrix;

    /**
     * @brief Evaluate the Post-Widder approximation of a continuous relaxation spectrum.
     *
     * For a relaxation function
     * \f[
     *   \Psi(t) = \int_0^\infty H(\tau)e^{-t/\tau}\,d\ln\tau,
     * \f]
     * this evaluates
     * \f[
     *   H_k(\tau) =
     *   \frac{(-k\tau)^k}{(k-1)!}\Psi^{(k)}(k\tau).
     * \f]
     *
     * @tparam k Post-Widder approximation order.
     * @param[in] psi Relaxation function whose spectrum is approximated.
     * @param[in] tau Relaxation time at which the spectrum is evaluated.
     * @return Approximate relaxation-spectrum value.
     */
    template < int k >
    double evaluatePostWidderFormula( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > psi,
                                      double                                                                      tau )
    {
      autodiff::Real< k, double > evaluationTime( tau * k );
      const double                coefficient = pow( -tau * k, k ) / double( Marmot::Math::factorial( k - 1 ) );
      return coefficient * autodiff::derivatives( psi, autodiff::along( 1. ), autodiff::at( evaluationTime ) )[k];
    }

    /**
     * @brief Compute generalized-Maxwell branch moduli from a relaxation function.
     *
     * The continuous spectrum is discretized on logarithmically spaced relaxation
     * times using \f$E_\mu = \ln(s)H_k(\tau_\mu)\f$, where \f$s\f$ is the spacing.
     *
     * @tparam k Post-Widder approximation order.
     * @param[in] psi Relaxation function whose spectrum is approximated.
     * @param[in] relaxationTimes Branch relaxation times.
     * @param[in] spacing Ratio between adjacent relaxation times.
     * @param[in] gaussQuadrature Use a two-point Gauss rule within each logarithmic interval.
     * @return Elastic modulus of every Maxwell branch.
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
          elasticModuli( i ) = log( spacing ) * evaluatePostWidderFormula< k >( psi, tau );
        }
        else {
          elasticModuli( i ) = log( spacing ) / 2. *
                               ( evaluatePostWidderFormula< k >( psi, tau * pow( spacing, -sqrt( 3. ) / 6. ) ) +
                                 evaluatePostWidderFormula< k >( psi, tau * pow( spacing, sqrt( 3. ) / 6. ) ) );
        }
      }
      return elasticModuli;
    }

    /**
     * @brief Generate logarithmically spaced Maxwell-branch relaxation times.
     * @param n Number of Maxwell branches.
     * @param min First relaxation time.
     * @param spacing Ratio between adjacent relaxation times.
     * @return Relaxation times of the Maxwell branches.
     */
    Properties generateRelaxationTimes( int n, double min, double spacing );

    /**
     * @brief Update branch state variables for one incremental strain step.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars In/out branch state matrix (6 x nMaxwell).
     * @param dStrain Strain increment in Voigt notation.
     * @param unitD_ijkl Elastic stiffness matrix used to map increment to stress-like branch updates.
     */
    void updateStateVarMatrix( const double                 dT,
                               const Properties&            elasticModuli,
                               const Properties&            relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStrain,
                               const Marmot::Matrix6d&      unitD_ijkl );

    /**
     * @brief Evaluate accumulated incremental Wiechert response from branch state variables.
     *
     * @details
     * Adds branch contributions to `uniaxialStiffness` and `dStress` (outputs are incremented, not reset).
     *
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars Branch state matrix (6 x nMaxwell).
     * @param uniaxialStiffness In/out accumulated scalar stiffness contribution.
     * @param dStress In/out accumulated stress increment in Voigt notation.
     * @param factor Optional scaling factor for all branch contributions.
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
