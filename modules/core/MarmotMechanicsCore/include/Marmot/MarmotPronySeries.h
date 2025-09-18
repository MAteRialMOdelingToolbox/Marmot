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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::ContinuumMechanics::Viscoelasticity {

  namespace PronySeries {
    using namespace Marmot;
    using namespace Eigen;
/**
 * @brief Material properties for a generalized Maxwell (Prony series) model.
 *
 * This structure stores the parameters required for evaluating the viscoelastic
 * stress response of a material represented by a Prony series expansion.
 */
    struct Properties {
   /**
   * @brief Number of Prony series terms.
   *
   * Corresponds to the number of Maxwell branches in the generalized Maxwell model.
   */
      size_t                  nPronyTerms;
  /**
   * @brief Ultimate (equilibrium) stiffness matrix.
   *
   * Represents the stiffness contribution that persists at long times
   * (\f$t \to \infty\f$), outside the viscoelastic relaxation terms.
   */
      Matrix6d                ultimateStiffnessMatrix;
  /**
   * @brief Stiffness matrices of each Prony series term.
   *
   * Stored as a block matrix of size \f$6 \times (6 \cdot nPronyTerms)\f$,
   * with each \f$6\times 6\f$ block corresponding to the stiffness matrix
   * of one Maxwell branch.
   */
      Matrix< double, 6, -1 > pronyStiffnesses;
  /**
   * @brief Relaxation times of each Prony series term.
   *
   * Stored as a block matrix of size \f$6 \times (6 \cdot nPronyTerms)\f$,
   * with each \f$6\times 6\f$ block containing the relaxation times
   * associated with one Maxwell branch.
   */
      Matrix< double, 6, -1 > pronyRelaxationTimes;
    };
/**
 * @typedef StateVarMatrix
 * @brief Matrix of viscoelastic state variables.
 *
 * Dynamic-size matrix of shape \f$[6 \times n]\f$, where each
 * \f$6\times 1\f$ block corresponds to the internal state variables
 * of a Maxwell branch in the Prony series.
 */
    typedef Eigen::Matrix< double, 6, Eigen::Dynamic > StateVarMatrix;
/**
 * @typedef mapStateVarMatrix
 * @brief Mapped view of a state variable matrix.
 *
 * Provides an Eigen::Map wrapper around a contiguous block of memory
 * that stores state variables in the same layout as StateVarMatrix,
 * avoiding copies.
 */
    typedef Eigen::Map< StateVarMatrix >               mapStateVarMatrix;
/**
 * @brief Evaluate the Prony series viscoelastic response.
 *
 * Computes the stress and tangent stiffness matrix contributions of a
 * generalized Maxwell (Prony series) viscoelastic material, and optionally
 * updates the internal state variables.
 *
 * The stress is updated as
 * \f[
 *   \sigma^{k+1}_{ij} \;+=\; C^\infty_{ijkl} \, \Delta \varepsilon_{kl}
 *          \;+\; \sum_{\mu=1}^{M} \left[
 *              \frac{\eta^\mu - \eta^\mu e^{-\Delta t / \tau^\mu}}{\Delta t}\,\Delta\varepsilon_{kl}
 *              \;-\; (\sigma^{ve,k}_{kl} - e^{-\Delta t/\tau^\mu} \sigma^{ve,k}_{kl})
 *          \right],
 * \f]
 * where \f$C^\infty_{ijkl}\f$ is the ultimate stiffness, \f$\tau^\mu\f$ the relaxation
 * times, \f$\eta^\mu = \tau^\mu C^\mu_{ijkl}\f$ the effective viscosities, and \f$\sigma^{ve,k}_{kl}\f$
 * the internal state variables.
 *
 * @param[in] props Prony series material properties (relaxation times,
 *                  stiffnesses, number of terms, and ultimate stiffness).
 * @param[in,out] stress 6-component stress vector to be updated.
 * @param[in,out] stiffness 6×6 tangent stiffness matrix to be updated.
 * @param[in,out] stateVars State variable matrix of size [6×(6·nPronyTerms)],
 *                          storing internal viscoelastic strain-like variables.
 * @param[in] dStrain Incremental strain vector (6-components).
 * @param[in] dT Time increment.
 * @param[in] updateStateVars If true, update the state variables for the next step.
 */
    void evaluatePronySeries( const Properties&               props,
                              Vector6d&                       stress,
                              Matrix6d&                       stiffness,
                              Eigen::Ref< mapStateVarMatrix > stateVars,
                              const Vector6d&                 dStrain,
                              const double                    dT,
                              const bool                      updateStateVars = false );
/**
 * @brief Update the Prony series state variables.
 *
 * Propagates the internal viscoelastic state variables
 * \f$\sigma^{ve,k}_{kl}\f$ from time step \f$t^k\f$ to \f$t^{k+1}\f$ according to
 * \f[
 *   \sigma^{ve,k+1}_{kl} = e^{-\Delta t / \tau^\mu}\, \sigma^{ve,k}_{kl}
 *             + \left( \frac{\eta^\mu}{\Delta t}(1 - e^{-\Delta t / \tau^\mu}) \right)\Delta\varepsilon_{kl},
 * \f]
 * where \f$\tau^\mu\f$ are the relaxation times and \f$\eta^\mu = \tau^\mu C^\mu_{ijkl}\f$.
 *
 * @param[in] props Prony series material properties (relaxation times,
 *                  stiffnesses, and number of terms).
 * @param[in,out] stateVars State variable matrix of size [6×(6·nPronyTerms)],
 *                          storing internal viscoelastic strain-like variables.
 * @param[in] dStrain Incremental strain vector (6-components).
 * @param[in] dT Time increment.
 */

    void updateStateVars( const Properties&               props,
                          Eigen::Ref< mapStateVarMatrix > stateVars,
                          const Vector6d&                 dStrain,
                          const double                    dT );

  } // namespace PronySeries
} // namespace Marmot::ContinuumMechanics::Viscoelasticity
