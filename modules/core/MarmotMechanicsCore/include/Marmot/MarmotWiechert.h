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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @namespace Wiechert
   * @brief Utilities for branch-wise Wiechert/Kelvin-chain viscoelastic updates in 3D Voigt notation.
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

    // TODO(v26.05): implement optional Post-Widder-based branch modulus generation helpers
    // (e.g. computeElasticModuli_* / generateRelaxationTimes) when calibration workflow needs them.
    // Right now only a Maxwell Element is implemented.

    /**
     * @brief Initialize branch-wise elastic moduli with a constant value.
     * @param nMaxwell Number of Maxwell branches.
     * @param n Elastic modulus assigned to every branch.
     * @return Vector of size `nMaxwell` with all entries equal to `n`.
     */
    Properties initializeElasticModuli( int nMaxwell, double n );

    /**
     * @brief Initialize branch-wise relaxation times with a constant value.
     * @param nMaxwell Number of Maxwell branches.
     * @param m Relaxation time assigned to every branch.
     * @return Vector of size `nMaxwell` with all entries equal to `m`.
     */
    Properties initializeRelaxationTimes( int nMaxwell, double m );

    /**
     * @brief Update branch state variables for one incremental strain step.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars In/out branch state matrix (6 x nMaxwell).
     * @param dStress Incremental driving quantity in Voigt notation.
     * @param unitD_ijkl Elastic stiffness matrix used to map increment to stress-like branch updates.
     */
    void updateStateVarMatrix( const double                 dT,
                               const Properties&            elasticModuli,
                               const Properties&            relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStress,
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
