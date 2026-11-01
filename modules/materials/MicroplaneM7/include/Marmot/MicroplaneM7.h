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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MicroplaneM7Formulation.h"

namespace Marmot::Materials {
  /**
   * @brief Implementation of the microplane model M7 for 3D stress states.
   *
   * @details The base formulation follows Caner & Bažant (2013), "Microplane model M7 for
   * plain concrete. I: Formulation" and "II: Calibration and verification", Journal of
   * Engineering Mechanics, 139(12). The closed-form transition function of the elastic unloading
   * stiffness (see MicroplaneM7Formulation::transitionFunctionClosedForm), improving convergence behavior
   * in implicit simulations, is introduced in Basche et al. (2026). The constitutive laws are implemented in
   * #MicroplaneM7Formulation, which is shared with the gradient-enhanced extension #GM7.
   */
  class MicroplaneM7 : public MarmotMaterialHypoElastic {
  public:
    MicroplaneM7( const double* materialProperties, int nMaterialProperties, int materialNumber );

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
      // fractions of the microplanes on the tensile, compressive and shear boundary, for postprocessing only
      stateLayout.add( "activeTensileBoundary", 1 );
      stateLayout.add( "activeCompressiveBoundary", 1 );
      stateLayout.add( "activeShearBoundary", 1 );
      stateLayout.add( "sig_N0", nqp );
      stateLayout.add( "sig_M0", nqp );
      stateLayout.add( "sig_L0", nqp );
      stateLayout.add( "eps_N0pos", nqp );
      stateLayout.add( "eps_N0neg", nqp );
      stateLayout.add( "strain_macro0", 6 );
      stateLayout.finalize();
    };

    double getDensity( const double* stateVars ) const override;

  protected:
    /// @brief The constitutive laws of M7, evaluated with the local strain.
    const MicroplaneM7Formulation m7_;

    using NormalResponse = MicroplaneM7Formulation::NormalResponse;

    /// @brief Per-microplane stress and macroscopically-projected tangent, stacked over all #nqp microplanes.
    struct allmp_response {
      nmp_sized  stress  = nmp_sized::Zero();
      nmp6_sized tangent = nmp6_sized::Zero();
    };

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    /// @brief Performs the full microplane stress update (kinematic projection, per-microplane
    /// normal/shear boundary evaluation, and homogenization back to macroscopic stress/tangent)
    /// for one strain increment, committing the resulting state variables into @p stateVars.
    /// @return {macroscopic stress, analytical tangent}.
    std::pair< Vector6d, Matrix6d > stress_update( const Vector6d& dE_macro, double* stateVars ) const;

    /// @brief Evaluates the normal-direction response of a single microplane: compares the
    /// elastic prediction (MicroplaneM7Formulation::normal_elasticLaw) against the tensile
    /// (MicroplaneM7Formulation::tensile_boundary) and compressive
    /// (MicroplaneM7Formulation::deviatoric_boundary + MicroplaneM7Formulation::volumetric_boundary)
    /// boundaries and returns whichever governs, updating the corresponding state variables and
    /// boundary counters in place.
    NormalResponse computeMicroplaneNormalResponse( const size_t          qp_i,
                                                    const double          sig_N0,
                                                    const double          deps_N,
                                                    const double          eps_N,
                                                    const double          sig_v0,
                                                    const double          zeta,
                                                    const double          eps_d,
                                                    const double          eps_v,
                                                    const NormalResponse& volumetric_law,
                                                    double&               eps_N0pos,
                                                    double&               eps_N0neg,
                                                    int&                  tensileB_flag,
                                                    int&                  compB_flag ) const;
  };
} // namespace Marmot::Materials
