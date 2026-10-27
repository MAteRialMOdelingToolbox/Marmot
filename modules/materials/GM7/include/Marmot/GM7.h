/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Konstantin Basche konstantin.basche@uibk.ac.at
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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MicroplaneM7Formulation.h"

namespace Marmot::Materials {
  /**
   * @brief Nonlocal gradient-enhanced extension of the microplane model M7 for 3D stress states.
   *
   * @details The base formulation follows Caner & Bažant (2013), "Microplane model M7 for
   * plain concrete. I: Formulation" and "II: Calibration and verification", Journal of
   * Engineering Mechanics, 139(12). The tensile, deviatoric, volumetric and shear boundaries are
   * evaluated on an "over-nonlocal" strain (see #stress_update) that blends the local strain with a
   * nonlocally-averaged strain via the weighting parameter #m, with the nonlocal interaction
   * parameter #c, for regularized (mesh-independent) softening. Both the closed-form transition
   * function of the elastic unloading stiffness (see MicroplaneM7Formulation::transitionFunctionClosedForm),
   * improving convergence behavior in implicit simulations, and the nonlocal gradient-enhancement
   * of the boundaries, are introduced in
   * Basche et al. (2026). The constitutive laws of M7 are shared with #MicroplaneM7 through
   * #MicroplaneM7Formulation.
   */
  class GM7 : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 6 > {
  public:
    GM7( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /// @brief Number of microplanes; see MicroplaneM7Formulation::nqp.
    constexpr static size_t nqp = MicroplaneM7Formulation::nqp;

    using nmp_sized  = MicroplaneM7Formulation::nmp_sized;
    using mnmp_sized = MicroplaneM7Formulation::mnmp_sized;
    using nmp6_sized = MicroplaneM7Formulation::nmp6_sized;

    /// @brief Declares the per-microplane and scalar state variables managed by #stateLayout.
    void initializeStateLayout()
    {
      stateLayout.add( "sig_v0", 1 );
      stateLayout.add( "zeta", 1 );
      stateLayout.add( "activeTensileBoundary", 1 );
      stateLayout.add( "activeCompressiveBoundary", 1 );
      stateLayout.add( "activeShearBoundary", 1 );
      stateLayout.add( "sig_N0", nqp );
      stateLayout.add( "sig_M0", nqp );
      stateLayout.add( "sig_L0", nqp );
      stateLayout.add( "eps_N0pos", nqp );
      stateLayout.add( "eps_N0neg", nqp );
      stateLayout.add( "strain_macro0OverNonLocal", 6 );
      stateLayout.add( "strain_macro0Local", 6 );
      stateLayout.finalize();
    };

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  protected:
    /// @brief The constitutive laws of M7; the elastic prediction is evaluated with the local, the
    /// boundaries with the over-nonlocal strain.
    const MicroplaneM7Formulation m7_;

    const double c; ///< @brief Nonlocal interaction parameter (length squared, e.g. mm²).
    const double m; ///< @brief Local/nonlocal weighting parameter for the "over-nonlocal" strain used to evaluate the
                    ///< boundaries (1.0 = purely nonlocal, 0.0 = purely local); see #stress_update.

    /// @brief Per-microplane normal-direction response with the tangent already projected to
    /// macroscopic (Voigt) space, split into the tangent w.r.t. the local strain increment
    /// (#tangent_l) and w.r.t. the over-nonlocal strain increment (#tangent_onl).
    struct NormalResponse {
      double   stress;
      Vector6d tangent_l;
      Vector6d tangent_onl;
    };

    /// @brief Per-microplane stress and macroscopically-projected local/over-nonlocal tangents, stacked over all #nqp
    /// microplanes.
    struct allmp_response {
      nmp_sized  stress      = nmp_sized::Zero();
      nmp6_sized tangent_l   = nmp6_sized::Zero();
      nmp6_sized tangent_onl = nmp6_sized::Zero();
    };

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    /// @brief Performs the full microplane stress update (kinematic projection, per-microplane
    /// normal/shear boundary evaluation, and homogenization back to macroscopic stress and
    /// tangents) for one strain increment, always committing the resulting per-microplane state
    /// variables into @p stateVars.
    /// @param dE_macroLocal Local macroscopic strain increment.
    /// @param dE_macroOverNonLocal Over-nonlocal macroscopic strain increment, `m * dE_macroNonLocal
    /// + (1 - m) * dE_macroLocal`; used to evaluate the tensile, deviatoric, volumetric and
    /// shear boundaries.
    /// @return {macroscopic stress, tangent w.r.t. local strain, tangent w.r.t. over-nonlocal strain}.
    std::tuple< Vector6d, Matrix6d, Matrix6d > stress_update( const Vector6d& dE_macroLocal,
                                                              const Vector6d& dE_macroOverNonLocal,
                                                              double*         stateVars ) const;

    /// @brief Evaluates the normal-direction response of a single microplane: compares the
    /// (local-strain) elastic prediction (MicroplaneM7Formulation::normal_elasticLaw) against the
    /// tensile (MicroplaneM7Formulation::tensile_boundary) and compressive
    /// (MicroplaneM7Formulation::deviatoric_boundary + MicroplaneM7Formulation::volumetric_boundary,
    /// evaluated on the over-nonlocal strain) boundaries and returns whichever governs, updating
    /// the corresponding state variables and boundary-exceeded counters in place.
    NormalResponse computeMicroplaneNormalResponse( const size_t                                   qp_i,
                                                    const double                                   sig_N0,
                                                    const double                                   deps_N_l,
                                                    const double                                   eps_N_l,
                                                    const double                                   eps_N_onl,
                                                    const double                                   eps_d_onl,
                                                    const double                                   eps_v_onl,
                                                    const double                                   sig_v0,
                                                    const double                                   zeta,
                                                    const MicroplaneM7Formulation::NormalResponse& volumetric_law,
                                                    double&                                        eps_N0pos,
                                                    double&                                        eps_N0neg,
                                                    int&                                           tensileB_flag,
                                                    int&                                           compB_flag ) const;
  };
} // namespace Marmot::Materials
