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
#include "Marmot/MarmotMicroplane.h"
#include "Marmot/MarmotTypedefs.h"
#include <tuple>
#include <utility>

namespace Marmot::Materials {
  /**
   * @brief The constitutive laws of the microplane model M7, shared by #MicroplaneM7 and its
   * gradient-enhanced extension #GM7, which both hold an instance.
   *
   * @details The base formulation follows Caner & Bažant (2013), "Microplane model M7 for
   * plain concrete. I: Formulation" and "II: Calibration and verification", Journal of
   * Engineering Mechanics, 139(12). The closed-form transition function of the unloading
   * stiffness (see #transitionFunctionClosedForm) is introduced in Basche et al. (2026).
   *
   * All functions are pure: they neither read nor write state variables. The normal laws and
   * boundaries return the stress and its derivative w.r.t. the macroscopic strain they are
   * evaluated with. The shear law splits its tangents by the strain they are derivatives
   * with respect to: suffix `_elastic` for the strain of the elastic predictions (normal and
   * shear), suffix `_boundary` for the strain the boundaries are evaluated with. For
   * #MicroplaneM7 both are the same strain, and the full tangent is the sum of the two parts,
   * as all expressions are linear in the tangents passed in. #GM7 evaluates the elastic
   * predictions with the local and the boundaries with the over-nonlocal strain, so the two
   * parts are its tangents w.r.t. the local and the over-nonlocal strain.
   */
  class MicroplaneM7Formulation {
  public:
    using QuadratureScheme = Marmot::Microplane::QuadratureScheme;

    /// @brief The quadrature scheme this compile-time instance is built with. Fixed at compile
    /// time - changing it requires editing this value and recompiling, there is no runtime
    /// selection between schemes.
    constexpr static QuadratureScheme scheme = QuadratureScheme::N37;

    /// @brief Number of microplanes (integration points on the unit hemisphere) for #scheme.
    constexpr static size_t nqp = Marmot::Microplane::nqpFromScheme( scheme );

    typedef Eigen::Matrix< double, 1, nqp >               nmp_sized;
    typedef Eigen::Map< Eigen::Matrix< double, 1, nqp > > mnmp_sized;
    typedef Eigen::Matrix< double, 6, nqp >               nmp6_sized;

    /// @brief The free (calibratable) material parameters.
    struct Parameters {
      double E;                  ///< @brief Young's modulus.
      double nu;                 ///< @brief Poisson's ratio.
      double k1;                 ///< @brief Radial scaling parameter (typical value 1.5e-4).
      double k2;                 ///< @brief Horizontal asymptote of the frictional boundary (typical value 110).
      double k3;                 ///< @brief Shape of the volumetric boundary (typical value 30).
      double k4;                 ///< @brief Shape of the volumetric boundary (typical value 100).
      double k5;                 ///< @brief Volumetric-deviatoric coupling at low pressures (typical value 1e-4).
      double c3;                 ///< @brief Slope of the postpeak in uniaxial tension (default 4.0).
      double fc;                 ///< @brief Compressive strength in MPa, for the deviatoric boundary (default 36.0).
      double transitionWidth;    ///< @brief Stress offset from the virtual elastic loading path where the
                                 ///< unloading transition starts (default 4.0).
      double transitionExponent; ///< @brief Exponent of the unloading transition function (default 0.4).
    };

    /// @brief Reads E, nu, k1 ... k5 from the material properties 0 to 6, and the optional c3, fc,
    /// transitionWidth and transitionExponent from the properties @p firstOptional to
    /// @p firstOptional + 3, using their default values for the ones not given. All properties in
    /// front of @p firstOptional are required.
    /// @throws std::invalid_argument if fewer than @p firstOptional material properties are given.
    static Parameters readParameters( const double* materialProperties, int nMaterialProperties, int firstOptional );

    explicit MicroplaneM7Formulation( const Parameters& parameters );

    /// @brief The free material parameters.
    const Parameters p;

    // Fixed material parameters of the M7 model according to Caner and Bazant (2013)
    // clang-format off
    const double fc0_reference = 15.08;   ///< @brief Reference compressive strength in MPa.
    const double E0_reference  = 20.0e3;  ///< @brief Reference Young's modulus in MPa.
    const double c1            = 8.9e-2;  ///< @brief Controls the uniaxial tensile strength.
    const double c2            = 17.6e-2; ///< @brief Controls the roundness of the peak in uniaxial tension.
    const double c4            = 50.0;    ///< @brief Controls the slope of the postpeak tail in uniaxial compression.
    const double c5            = 3500.0;  ///< @brief Controls the volumetric expansion under compression.
    const double c6            = 20.0;    ///< @brief Controls the roundness of the peak in volumetric expansion under compression.
    const double c7            = 1.0;     ///< @brief Controls the slope of the initial postpeak in uniaxial compression.
    const double c8            = 8.0;     ///< @brief Controls the peak strength in uniaxial compression.
    const double c9            = 1.2e-2;  ///< @brief Controls the peak roundness in uniaxial compression.
    const double c10           = 0.33;    ///< @brief Controls the effective friction coefficient.
    const double c11           = 0.5;     ///< @brief Initial cohesion in frictional response.
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
    // clang-format on

    /// @brief Kinematic (N-T split) projection tensors and integration weights for all #nqp
    /// microplanes; see Marmot::Microplane::NTSplit.
    const Marmot::Microplane::NTSplit< scheme > ntSplit;

    /// @brief Undamaged microplane normal stiffness E / (1 - 2 nu).
    double normalStiffness() const;

    /// @brief Undamaged microplane tangential stiffness; there is no damage of the shear stiffness.
    double shearStiffness() const;

    /// @brief Microplane normal stress and its derivative w.r.t. the macroscopic strain, i.e.
    /// `tangent = d(stress)/d(microplane strain) * projection.N.col(qp_i)`.
    struct NormalResponse {
      double   stress;
      Vector6d tangent;
    };

    /// @brief Elastic (undamaged) stress prediction for one shear direction; see #shear_elasticLaw.
    struct ElasticPrediction {
      double stress;
      double tangent;
    };

    /// @brief Quantities of the shear law up to the point where the R27/R28 frictional-boundary
    /// branches are combined into a single shear boundary.
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
      Vector6d          dsigDachddStrain_boundary;
      Vector6d          dR27ddStrain_elastic;
      Vector6d          dR27ddStrain_boundary;
      Vector6d          dR28ddStrain_boundary;
    };

    /// @brief Microplane shear response (both in-plane directions L and M), with the tangent of
    /// each direction split into the derivatives w.r.t. the strain of the elastic predictions
    /// (`_elastic`) and w.r.t. the strain the boundaries are evaluated with (`_boundary`).
    struct ShearResponse {
      double   stressL;
      double   stressM;
      Vector6d tangentL_elastic;
      Vector6d tangentL_boundary;
      Vector6d tangentM_elastic;
      Vector6d tangentM_boundary;
      bool     boundaryExceeded; ///< @brief True if the shear boundary governs.
    };

    /// @brief Normal-direction elastic/damage law (tension and compression branches), including
    /// the transition of the elastic unloading stiffness to the undamaged stiffness in tension
    /// via #transitionFunctionClosedForm.
    NormalResponse normal_elasticLaw( const double sig_N0,
                                      const double deps_N,
                                      const double eps_N,
                                      const double eps_N0pos,
                                      const double eps_N0neg,
                                      const double sig_v0,
                                      const double zeta,
                                      const size_t qp_i ) const;

    /// @brief Closed-form integration of the unloading transition from damaged to undamaged
    /// stiffness, used by #normal_elasticLaw.
    std::pair< double, double > transitionFunctionClosedForm( const double eps_N0,
                                                              const double deps_N,
                                                              const double sig_N0,
                                                              const double E_N,
                                                              const double E_N0 ) const;

    /// @brief Blends the damaged stiffness E_N and the undamaged stiffness E_N0 over the
    /// unloading transition, returning the regularized stiffness.
    double unloadingTransitionStiffness( const double s_gap, const double E_N, const double E_N0 ) const;

    /// @brief Volumetric compressive boundary (same for all microplanes) for the macroscopic strain @p strain_macro.
    NormalResponse volumetric_boundary( const double sig_v0, const Vector6d& strain_macro ) const;

    /// @brief Deviatoric compressive boundary for microplane @p qp_i.
    NormalResponse deviatoric_boundary( const double eps_d, const double eps_v, const size_t qp_i ) const;

    /// @brief Tensile boundary for microplane @p qp_i.
    NormalResponse tensile_boundary( const double sig_v0, const double eps_N, const size_t qp_i ) const;

    /// @brief Shear law according to Caner and Bazant (2013): the elastic shear prediction is
    /// limited by the frictional boundary, switching between its R27/R28 branches on the sign of
    /// the normal stress.
    /// @param eps_v Volumetric strain the shear boundary is evaluated with.
    /// @param sig_N Microplane normal stress.
    /// @param sig_N_tangent_elastic Derivative of @p sig_N w.r.t. the strain of the elastic prediction.
    /// @param sig_N_tangent_boundary Derivative of @p sig_N w.r.t. the strain the boundaries are evaluated with.
    ShearResponse shear_law( const double    eps_v,
                             const double    sig_N,
                             const double    sig_M0,
                             const double    sig_L0,
                             const double    deps_M,
                             const double    deps_L,
                             const Vector6d& sig_N_tangent_elastic,
                             const Vector6d& sig_N_tangent_boundary,
                             const size_t    qp_i ) const;

    /// @brief Elastic (undamaged) shear law for one in-plane direction (M or L).
    ElasticPrediction shear_elasticLaw( const double sig_S0, const double deps_S ) const;

  private:
    /// @brief Computes the quantities of #shear_law before the R27/R28 branches are combined.
    ShearBoundaryInputs computeShearBoundaryInputs( const double    eps_v,
                                                    const double    sig_N,
                                                    const double    sig_M0,
                                                    const double    sig_L0,
                                                    const double    deps_M,
                                                    const double    deps_L,
                                                    const Vector6d& sig_N_tangent_elastic,
                                                    const Vector6d& sig_N_tangent_boundary ) const;

    /// @brief Given the combined shear boundary and its tangents, returns the microplane shear
    /// stresses and tangents for both in-plane directions.
    ShearResponse computeShearResponseFromBoundary( const ShearBoundaryInputs& in,
                                                    const double               shear_boundary,
                                                    const Vector6d&            shearBoundaryTangent_elastic,
                                                    const Vector6d&            shearBoundaryTangent_boundary,
                                                    const size_t               qp_i ) const;
  };
} // namespace Marmot::Materials
