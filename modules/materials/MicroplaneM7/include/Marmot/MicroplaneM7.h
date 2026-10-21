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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMicroplane.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {
  /**
   * @brief Implementation of the microplane model M7 for 3D stress states.
   *
   * @details The base formulation follows Caner & Bažant (2013), "Microplane model M7 for
   * plain concrete. I: Formulation" and "II: Calibration and verification", Journal of
   * Engineering Mechanics, 139(12). The closed-form normal-boundary transition function
   * (see #transitionFunctionClosedForm) improving convergence behavior in implicit
   * simulations is introduced in Basche et al. (2026).
   */
  class MicroplaneM7 : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    MicroplaneM7( const double* materialProperties, int nMaterialProperties, int materialNumber );

    using QuadratureScheme = Marmot::Microplane::QuadratureScheme;

    /// @brief The quadrature scheme this compile-time instance is built with. Fixed at compile
    /// time - changing it requires editing this value and recompiling, there is no runtime
    /// selection between schemes.
    constexpr static QuadratureScheme scheme = QuadratureScheme::N37;

    /// @brief Number of microplanes (integration points on the unit hemisphere) for #scheme.
    constexpr static size_t nqp = Marmot::Microplane::nqpFromScheme( scheme );

    /// @brief If true, #computeStress uses the analytical tangent; otherwise a central-difference
    /// approximation is used instead.
    bool analyticalTangent = true;

    typedef Eigen::Matrix< double, 1, nqp >               nmp_sized;
    typedef Eigen::Map< Eigen::Matrix< double, 1, nqp > > mnmp_sized;
    typedef Eigen::Matrix< double, 6, nqp >               nmp6_sized;

    /// @brief Declares the per-microplane and scalar state variables managed by #stateLayout.
    void initializeStateLayout()
    {
      stateLayout.add( "sig_v0", 1 );
      stateLayout.add( "zeta", 1 );
      stateLayout.add( "sig_N0", nqp );
      stateLayout.add( "sig_M0", nqp );
      stateLayout.add( "sig_L0", nqp );
      stateLayout.add( "eps_N0pos", nqp );
      stateLayout.add( "eps_N0neg", nqp );
      // further state variables are added only for postprocessing, nor used in the stress update itself
      stateLayout.add( "shear_boundary", nqp );
      stateLayout.add( "normal_boundary", nqp );
      stateLayout.add( "strain_macro0", 6 );
      stateLayout.finalize();
    };

    double getDensity( const double* stateVars ) const override;

  protected:
    // Free (calibratable) material parameters
    const double& E;  ///< @brief Young's modulus.
    const double& nu; ///< @brief Poisson's ratio.
    const double& k1; ///< @brief Radial scaling parameter (typical value 1.5e-4).
    const double&
      k2; ///< @brief Controls the horizontal asymptote value in the frictional boundary (typical value 110).
    const double& k3; ///< @brief Controls the shape of the volumetric boundary (typical value 30).
    const double& k4; ///< @brief Controls the shape of the volumetric boundary (typical value 100).
    const double& k5; ///< @brief Controls the volumetric-deviatoric coupling at low pressures (typical value 1e-4).

    // Fixed material parameters of the M7 model according to Caner and Bazant (2013)
    // clang-format off
    const double fc0_reference = 15.08;   ///< @brief Reference compressive strength in MPa.
    const double E0_reference  = 20.0e3;  ///< @brief Reference Young's modulus in MPa.
    const double c1            = 8.9e-2;  ///< @brief Controls the uniaxial tensile strength.
    const double c2            = 17.6e-2; ///< @brief Controls the roundness of the peak in uniaxial tension (and the peak value).
    const double c3;                      ///< @brief Controls the slope of the postpeak in uniaxial tension; defaults to 4.0 (Caner and Bazant 2013) if not given as a material property.
    const double c4            = 50.0;    ///< @brief Controls the slope of the postpeak tail in uniaxial compression.
    const double c5            = 3500.0;  ///< @brief Controls the volumetric expansion under compression (no influence on the deviatoric boundary).
    const double c6            = 20.0;    ///< @brief Controls the roundness of the peak in volumetric expansion.
    const double c7            = 1.0;     ///< @brief Controls the slope of the initial postpeak in uniaxial compression.
    const double c8            = 8.0;     ///< @brief Controls the peak strength in uniaxial compression.
    const double c9            = 1.2e-2;  ///< @brief Controls the peak roundness in uniaxial compression.
    const double c10           = 0.33;    ///< @brief Controls the effective friction coefficient.
    const double c11           = 0.5;     ///< @brief Initial cohesion in the frictional response.
    const double c12           = 2.36;    ///< @brief Controls the change of cohesion with tensile volumetric strains.
    const double c13           = 4500.0;  ///< @brief Controls the unloading slope in tension.
    const double c14           = 300.0;   ///< @brief Controls the unloading slope at low hydrostatic compression.
    const double c15           = 4000.0;  ///< @brief Controls the transition from the unloading slope at high confinement to that at low confinement.
    const double c16           = 60.0;    ///< @brief Controls the unloading slope at high hydrostatic compression.
    const double c17           = 1.4;     ///< @brief Controls the tensile strength.
    const double c18           = 62.5;    ///< @brief [MPa] Controls the tensile cracking under compression; value per the correction in Caner and Bazant (2013).
    const double c19           = 1000.0;  ///< @brief Controls the tensile softening rate induced by compression.
    const double c20           = 1.8;     ///< @brief Controls the volumetric-deviatoric coupling at high pressures.
    const double c21           = 250.0;   ///< @brief [MPa]
    /// @brief Compressive strength in MPa, needed for the deviatoric boundary.
    const double fc;
    // clang-format on

    /// @brief Stress offset from the virtual elastic loading path where the unloading transition starts.
    const double transitionWidth;
    /// @brief Exponent controlling the shape of the unloading transition function; see #transitionFunctionClosedForm.
    const double transitionExponent;

    /// @brief Per-microplane normal-direction response with the tangent already projected to
    /// macroscopic (Voigt) space, i.e. `tangent = d(stress)/d(microplane strain) * projection.N.col(qp_i)`.
    /// Used by #normal_elasticLaw, #tensile_boundary, #deviatoric_boundary and #volumetric_boundary.
    struct NormalResponse {
      double   stress;
      Vector6d tangent;
    };

    /// @brief Per-microplane stress and macroscopically-projected tangent, stacked over all #nqp microplanes.
    struct allmp_response {
      nmp_sized  stress  = nmp_sized::Zero();
      nmp6_sized tangent = nmp6_sized::Zero();
    };

    /// @brief Elastic (undamaged) stress prediction for one shear direction; see #shear_elasticLaw.
    struct ElasticPrediction {
      double stress;
      double tangent;
    };

    /// @brief Quantities shared by #shear_lawNguyen and #shear_lawCaner up to the point where
    /// they diverge on how they combine the R27/R28 frictional-boundary branches.
    struct ShearBoundaryInputs {
      double            E_T;
      double            sig_dach;
      ElasticPrediction elasticL_prediction;
      ElasticPrediction elasticM_prediction;
      double            shear_stress_elastic;
      double            R27_denom;
      double            R27;
      double            R28_denom;
      double            R28;
      Vector6d          dsigDachddStrain;
      Vector6d          dR27ddStrain;
      Vector6d          dR28ddStrain;
    };

    /// @brief Kinematic (N-T split) projection tensors and integration weights for all #nqp
    /// microplanes; see Marmot::Microplane::NTSplit.
    Marmot::Microplane::NTSplit< scheme > ntSplit_;

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    /// @brief Performs the full microplane stress update (kinematic projection, per-microplane
    /// normal/shear boundary evaluation, and homogenization back to macroscopic stress/tangent)
    /// for one strain increment.
    /// @tparam T Strain-increment type; accepts both `Vector6d` and `Eigen::Map<const Vector6d>`.
    /// @brief Performs one stress update on the given @p stateVars buffer, always committing the
    /// resulting per-microplane state variables into it. #computeStress passes either the real
    /// state (for the single analytical-tangent evaluation) or a disposable scratch copy of it
    /// (for the central-difference probing evaluations used when #analyticalTangent is false),
    /// so this function itself never needs to know which case it is.
    template < typename T >
    std::pair< Vector6d, Matrix6d > stress_update( const T& dE_macro, double* stateVars ) const;

    /// @brief Evaluates the normal-direction response of a single microplane: compares the
    /// elastic prediction (#normal_elasticLaw) against the tensile (#tensile_boundary) and
    /// compressive (#deviatoric_boundary + #volumetric_boundary) boundaries and returns
    /// whichever governs, updating the corresponding state variables in place.
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
                                                    double&               normal_boundary ) const;

    /// @brief Normal-direction elastic/damage law (tension and unloading branches), including
    /// the transition to the tensile boundary via #transitionFunctionClosedForm.
    NormalResponse normal_elasticLaw( const double sig_N0,
                                      const double deps_N,
                                      const double eps_N,
                                      const double eps_N0pos,
                                      const double eps_N0neg,
                                      const double sig_v0,
                                      const double zeta,
                                      const size_t qp_i ) const;

    /// @brief Undamaged normal-direction elastic law. Kept as a documented alternative to
    /// #normal_elasticLaw (superseded); not called from #stress_update.
    NormalResponse normal_elasticLawUndamaged( const double sig_N0, const double deps_N, const size_t qp_i ) const;

    /// @brief Sub-stepped (explicit) integration of the unloading transition from damaged to undamaged stiffness.
    /// Kept as a documented alternative to #transitionFunctionClosedForm (superseded, replaced
    /// for better computational efficiency); not called from #stress_update.
    std::pair< double, double > transitionFunction( const double eps_N0,
                                                    const double deps_N,
                                                    const double sig_N0,
                                                    const double E_N,
                                                    const double E_N0 ) const;

    /// @brief Closed-form integration of the unloading transition from damaged to undamaged stiffness,
    /// active counterpart to #transitionFunction, used by #normal_elasticLaw.
    std::pair< double, double > transitionFunctionClosedForm( const double eps_N0,
                                                              const double deps_N,
                                                              const double sig_N0,
                                                              const double E_N,
                                                              const double E_N0 ) const;

    /// @brief Blends the undamaged stiffness E_N and the damaged stiffness E_N0 over the
    /// unloading transition, returning the regularized stiffness and its derivatives w.r.t.
    /// strain and stress.
    std::tuple< double, double, double > unloadingTransitionStiffness( const double s_gap,
                                                                       const double E_N,
                                                                       const double E_N0 ) const;

    /// @brief Volumetric compressive boundary (same for all microplanes).
    NormalResponse volumetric_boundary( const double sig_v0, const Vector6d& strain_macro ) const;

    /// @brief Deviatoric compressive boundary for microplane @p qp_i.
    NormalResponse deviatoric_boundary( const double eps_d, const double eps_v, const size_t qp_i ) const;

    /// @brief Tensile boundary for microplane @p qp_i.
    NormalResponse tensile_boundary( const double sig_v0, const double eps_N, const size_t qp_i ) const;

    /// @brief Computes the quantities shared by #shear_lawNguyen and #shear_lawCaner (elastic
    /// shear prediction, the R27/R28 frictional-boundary branch values, and their strain
    /// derivatives) before they diverge on how R27/R28 are combined into a single shear boundary.
    ShearBoundaryInputs computeShearBoundaryInputs( const double    eps_v,
                                                    const double    sig_N,
                                                    const double    sig_M0,
                                                    const double    sig_L0,
                                                    const double    deps_M,
                                                    const double    deps_L,
                                                    const Vector6d& sig_N_tangent ) const;

    /// @brief Given the combined shear boundary and its tangent (computed by #shear_lawNguyen or
    /// #shear_lawCaner from a #ShearBoundaryInputs), returns the microplane shear stresses/tangents
    /// and updates the "shear_boundary" state variable.
    std::tuple< double, double, Vector6d, Vector6d > computeShearResponseFromBoundary( const ShearBoundaryInputs& in,
                                                                                       const double    shear_boundary,
                                                                                       const Vector6d& boundary_tangent,
                                                                                       const size_t    qp_i,
                                                                                       double* stateVars ) const;

    /// @brief Shear boundary with a smooth (tanh-blended) transition between the R27/R28
    /// frictional-boundary branches as proposed by Nguyen et al (2021). Kept as a documented alternative to
    /// #shear_lawCaner (superseded, unused in #stress_update).
    std::tuple< double, double, Vector6d, Vector6d > shear_lawNguyen( const double    eps_v,
                                                                      const double    sig_N,
                                                                      const double    sig_M0,
                                                                      const double    sig_L0,
                                                                      const double    deps_M,
                                                                      const double    deps_L,
                                                                      const Vector6d& sig_N_tangent,
                                                                      const size_t    qp_i,
                                                                      double*         stateVars ) const;

    /// @brief Shear boundary with a hard switch between the R27/R28 frictional-boundary
    /// branches on the sign of the normal stress; the active counterpart to #shear_lawNguyen,
    /// used by #stress_update.
    std::tuple< double, double, Vector6d, Vector6d > shear_lawCaner( const double    eps_v,
                                                                     const double    sig_N,
                                                                     const double    sig_M0,
                                                                     const double    sig_L0,
                                                                     const double    deps_M,
                                                                     const double    deps_L,
                                                                     const Vector6d& sig_N_tangent,
                                                                     const size_t    qp_i,
                                                                     double*         stateVars ) const;

    /// @brief Elastic (undamaged) shear law for one in-plane direction (M or L).
    ElasticPrediction shear_elasticLaw( const double sig_S0, const double deps_S ) const;
  };
} // namespace Marmot::Materials
