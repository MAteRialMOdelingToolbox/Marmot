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
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <fstream>
#include <functional>
#include <iostream>
//
namespace Marmot::Materials {

  /**
   * @namespace WiechertInterface
   * @brief Utility routines for interface-specific Wiechert (Kelvin-chain) viscoelastic updates.
   *
   * @details
   * The functions in this namespace provide:
   * - initialization helpers for Maxwell-branch properties,
   * - incremental evaluation of interface traction/stress contributions,
   * - state variable updates for each branch,
   * - robust computation of integration factors \f$\lambda\f$ and \f$\beta\f$.
   */
  namespace WiechertInterface {

    /// @brief Dynamic vector of branch-wise material properties (e.g. moduli or relaxation times).
    typedef Eigen::VectorXd Properties;

    /// @brief Non-owning map view onto a `Properties` buffer.
    typedef Eigen::Map< Properties > mapProperties;

    /// @brief State variable matrix for \f$\Delta \mathbf{f}_{uu}\f$-type branch forces (3 x nMaxwell).
    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrix_force_uu;

    /// @brief State variable matrix for \f$\Delta \mathbf{f}_{u\bar{\varepsilon}}\f$-type branch forces (3 x nMaxwell).
    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrix_force_us;

    /// @brief State variable matrix for \f$\Delta \boldsymbol{\sigma}^{\mathrm{surf}}_{Z}\f$ terms (9 x nMaxwell).
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_Z;

    /// @brief State variable matrix for \f$\Delta \boldsymbol{\sigma}^{\mathrm{surf}}_{Y}\f$ terms (9 x nMaxwell).
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_Y;

    /// @brief State variable matrix for coupling surface stresses \f$\Delta \boldsymbol{\sigma}^{\mathrm{surf}}_{us}\f$
    /// (9 x nMaxwell).
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_surface_stress_us;

    /// @brief Non-owning map view of `StateVarMatrix_force_uu`.
    typedef Eigen::Map< StateVarMatrix_force_uu > mapStateVarMatrix_force_uu;

    /// @brief Non-owning map view of `StateVarMatrix_force_us`.
    typedef Eigen::Map< StateVarMatrix_force_us > mapStateVarMatrix_force_us;

    /// @brief Non-owning map view of `StateVarMatrix_surface_stress_Z`.
    typedef Eigen::Map< StateVarMatrix_surface_stress_Z > mapStateVarMatrix_surface_stress_Z;

    /// @brief Non-owning map view of `StateVarMatrix_surface_stress_Y`.
    typedef Eigen::Map< StateVarMatrix_surface_stress_Y > mapStateVarMatrix_surface_stress_Y;

    /// @brief Non-owning map view of `StateVarMatrix_surface_stress_us`.
    typedef Eigen::Map< StateVarMatrix_surface_stress_us > mapStateVarMatrix_surface_stress_us;

    // TODO(v26.05): Support generalized Maxwell chains with distinct branch moduli and relaxation times,
    // supplied directly or generated from an interface-specific approximation such as a power law.

    /**
     * @brief Update branch state variables for the `force_uu` contribution.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_force_uu In-place state variable matrix (3 x nMaxwell).
     * @param dforce_uu Incremental jump-like quantity driving this branch contribution.
     * @param unitH_inv_ij Geometric/material coupling matrix.
     */
    void updateStateVarMatrix_force_uu( const double                          dT,
                                        const Properties&                     elasticModuli,
                                        const Properties&                     relaxationTimes,
                                        Eigen::Ref< StateVarMatrix_force_uu > stateVars_force_uu,
                                        const Marmot::Vector3d&               dforce_uu,
                                        const Marmot::Matrix3d&               unitH_inv_ij );

    /**
     * @brief Update branch state variables for the `force_us` contribution.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_force_us In-place state variable matrix (3 x nMaxwell).
     * @param dforce_us Incremental strain-like driving quantity.
     * @param unitH_inv_nF_ijk Coupling operator mapping to branch force space.
     */
    void updateStateVarMatrix_force_us( const double                          dT,
                                        const Properties&                     elasticModuli,
                                        const Properties&                     relaxationTimes,
                                        Eigen::Ref< StateVarMatrix_force_us > stateVars_force_us,
                                        const Marmot::Vector9d&               dforce_us,
                                        const Eigen::Matrix< double, 3, 9 >&  unitH_inv_nF_ijk );

    /**
     * @brief Update branch state variables for the surface-stress `Z` contribution.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_surface_stress_Z In-place state matrix (9 x nMaxwell).
     * @param dsurfaceStress_Z Incremental driving quantity.
     * @param unitZ_ijkl Surface projection/coupling tensor in matrix form.
     */
    void updateStateVarMatrix_surface_stress_Z(
      const double                                  dT,
      const Properties&                             elasticModuli,
      const Properties&                             relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_Z > stateVars_surface_stress_Z,
      const Marmot::Vector9d&                       dsurfaceStress_Z,
      const Marmot::Matrix9d&                       unitZ_ijkl );

    /**
     * @brief Update branch state variables for the surface-stress `Y` contribution.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_surface_stress_Y In-place state matrix (9 x nMaxwell).
     * @param dsurfaceStress_Y Incremental driving quantity.
     * @param unitYn_H_inv_Fn_ijkl Coupling tensor in matrix form.
     */
    void updateStateVarMatrix_surface_stress_Y(
      const double                                  dT,
      const Properties&                             elasticModuli,
      const Properties&                             relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_Y > stateVars_surface_stress_Y,
      const Marmot::Vector9d&                       dsurfaceStress_Y,
      const Marmot::Matrix9d&                       unitYn_H_inv_Fn_ijkl );

    /**
     * @brief Update branch state variables for coupled `surface_stress_us` terms.
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_surface_stress_us In-place state matrix (9 x nMaxwell).
     * @param djumpU Incremental displacement jump vector.
     * @param unitH_inv_nF_ijk Coupling operator mapping jump increments to stress space.
     */
    void updateStateVarMatrix_surface_stress_us(
      const double                                   dT,
      const Properties&                              elasticModuli,
      const Properties&                              relaxationTimes,
      Eigen::Ref< StateVarMatrix_surface_stress_us > stateVars_surface_stress_us,
      const Marmot::Vector3d&                        djumpU,
      const Eigen::Matrix< double, 3, 9 >&           unitH_inv_nF_ijk );

    /**
     * @brief Evaluate incremental Wiechert response contributions for interface quantities.
     *
     * @details
     * Accumulates branch-wise contributions into stiffness, force and surface-stress increments.
     * The output references are incremented (not reset), enabling additive assembly.
     *
     * @param dT Time increment.
     * @param elasticModuli Branch-wise elastic moduli, passed by const reference.
     * @param relaxationTimes Branch-wise relaxation times, passed by const reference.
     * @param stateVars_force_uu Branch state matrix for `force_uu` contribution.
     * @param stateVars_force_us Branch state matrix for `force_us` contribution.
     * @param stateVars_surface_stress_Z Branch state matrix for surface-stress `Z` terms.
     * @param stateVars_surface_stress_Y Branch state matrix for surface-stress `Y` terms.
     * @param stateVars_surface_stress_us Branch state matrix for coupled surface-stress terms.
     * @param uniaxialStiffness In/out accumulated uniaxial stiffness contribution.
     * @param dforce_uu In/out accumulated `force_uu` increment.
     * @param dforce_us In/out accumulated `force_us` increment.
     * @param dsurfaceStress_Z In/out accumulated surface-stress `Z` increment.
     * @param dsurfaceStress_Y In/out accumulated surface-stress `Y` increment.
     * @param dsurfaceStress_us In/out accumulated coupled surface-stress increment.
     * @param factor Optional scaling factor applied to all branch contributions.
     */
    void evaluateWiechert( const double                            dT,
                           const Properties&                       elasticModuli,
                           const Properties&                       relaxationTimes,
                           const StateVarMatrix_force_uu&          stateVars_force_uu,
                           const StateVarMatrix_force_us&          stateVars_force_us,
                           const StateVarMatrix_surface_stress_Z&  stateVars_surface_stress_Z,
                           const StateVarMatrix_surface_stress_Y&  stateVars_surface_stress_Y,
                           const StateVarMatrix_surface_stress_us& stateVars_surface_stress_us,
                           double&                                 uniaxialStiffness,
                           Marmot::Vector3d&                       dforce_uu,
                           Marmot::Vector3d&                       dforce_us,
                           Marmot::Vector9d&                       dsurfaceStress_Z,
                           Marmot::Vector9d&                       dsurfaceStress_Y,
                           Marmot::Vector9d&                       dsurfaceStress_us,
                           const double                            factor );

  } // namespace WiechertInterface
} // namespace Marmot::Materials
