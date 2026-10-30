#include "Marmot/MicroplaneM7Formulation.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::VoigtNotation;

  namespace {
    // Fallback values for the optional material properties (c3, fc, transitionWidth,
    // transitionExponent), used when they are not supplied by the caller.
    constexpr double defaultPostpeakTensionSlope = 4.0;  // c3, per Caner and Bazant (2013)
    constexpr double defaultCompressiveStrength  = 36.0; // fc [MPa], per Caner and Bazant (2013)
    constexpr double defaultTransitionWidth      = 4.0;
    constexpr double defaultTransitionExponent   = 0.4;

    // Safety floor on the normal-direction tangent, as a fraction of the undamaged stiffness
    // E_N0, to avoid a vanishing (or sign-flipping) tangent in normal_elasticLaw.
    constexpr double minStiffnessFraction = 0.02;

    // Tolerance used in normal_elasticLaw to detect that the tangent has already dropped
    // below the undamaged stiffness when deciding whether to (re-)enter the unloading-transition branch.
    constexpr double tangentTransitionTol = 1e-6;
  } // namespace

  MicroplaneM7Formulation::Parameters MicroplaneM7Formulation::readParameters( const double* materialProperties,
                                                                               int           nMaterialProperties,
                                                                               int           firstOptional )
  {
    if ( nMaterialProperties < firstOptional )
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": " << firstOptional << " material properties required, "
                                   << nMaterialProperties << " given." );

    const auto optional = [&]( int index, double defaultValue ) {
      return nMaterialProperties > index ? materialProperties[index] : defaultValue;
    };

    return { materialProperties[0],
             materialProperties[1],
             materialProperties[2],
             materialProperties[3],
             materialProperties[4],
             materialProperties[5],
             materialProperties[6],
             optional( firstOptional, defaultPostpeakTensionSlope ),
             optional( firstOptional + 1, defaultCompressiveStrength ),
             optional( firstOptional + 2, defaultTransitionWidth ),
             optional( firstOptional + 3, defaultTransitionExponent ) };
  }

  MicroplaneM7Formulation::MicroplaneM7Formulation( const Parameters& parameters ) : p( parameters ) {}

  const Marmot::Microplane::NTSplit< MicroplaneM7Formulation::scheme >& MicroplaneM7Formulation::ntSplit()
  {
    static const Marmot::Microplane::NTSplit< scheme > ntSplit_;
    return ntSplit_;
  }

  double MicroplaneM7Formulation::normalStiffness() const
  {
    return p.E / ( 1.0 - 2.0 * p.nu );
  }

  double MicroplaneM7Formulation::shearStiffness() const
  {
    return p.E / ( 1.0 - 2.0 * p.nu ) * ( 1.0 - 4.0 * p.nu ) / ( 1.0 + p.nu );
  }

  MicroplaneM7Formulation::NormalResponse MicroplaneM7Formulation::normal_elasticLaw( const double sig_N0,
                                                                                      const double deps_N,
                                                                                      const double eps_N,
                                                                                      const double eps_N0pos,
                                                                                      const double eps_N0neg,
                                                                                      const double sig_v0,
                                                                                      const double zeta,
                                                                                      const size_t qp_i ) const
  {
    const double E_N0 = normalStiffness(); // microplane normal undamaged stiffness
    double       stress;
    double       tangent;
    const double minStiffness = minStiffnessFraction * E_N0;

    if ( sig_N0 >= 0.0 ) {
      double fzeta       = 1.0 / ( 1.0 + 0.1 * zeta * zeta );
      tangent            = std::max( E_N0 * Math::exp( -c13 * eps_N0pos ) * fzeta, minStiffness );
      double sig_N       = sig_N0 + tangent * deps_N;
      double eps_N0      = eps_N - deps_N;
      double s_gap_trial = E_N0 * eps_N - sig_N;
      if ( s_gap_trial < p.transitionWidth && ( sig_N0 * deps_N ) < 0.0 &&
           std::abs( E_N0 - tangent ) > tangentTransitionTol ) {
        std::tie( stress, tangent ) = transitionFunctionClosedForm( eps_N0, deps_N, sig_N0, tangent, E_N0 );
        return { stress, tangent * ntSplit().pTensors.N.col( qp_i ) };
      }
      // ------------ original M7 (Caner 2013) --------------
      // if ( sig_N0 > ( E_N0 * eps_N ) && ( sig_N0 * deps_N ) < 0.0 ) {
      //   tangent = E_N0;
      // }
    }
    else {
      tangent = std::max( E_N0 * ( Math::exp( -c14 * std::abs( eps_N0neg ) /
                                              ( 1.0 + c15 * std::max( -sig_v0 / E_N0, 0.0 ) ) ) +
                                   c16 * std::max( -sig_v0 / E_N0, 0.0 ) ),
                          minStiffness );
    }

    stress = sig_N0 + tangent * deps_N;

    return { stress, tangent * ntSplit().pTensors.N.col( qp_i ) };
  }

  std::pair< double, double > MicroplaneM7Formulation::transitionFunctionClosedForm( const double eps_N0,
                                                                                     const double deps_N,
                                                                                     const double sig_N0,
                                                                                     const double E_N,
                                                                                     const double E_N0 ) const
  {
    const double transitionWidth    = p.transitionWidth;
    const double transitionExponent = p.transitionExponent;

    const double s_gap0 = E_N0 * eps_N0 - sig_N0;
    double       eps_N_start, deps_N_transition, s_gap_start;
    double       sig_N;
    if ( s_gap0 > transitionWidth ) { // transition not yet started -> calculate starting point
      eps_N_start              = ( transitionWidth - E_N * eps_N0 + sig_N0 ) / ( E_N0 - E_N );
      deps_N_transition        = deps_N - ( eps_N_start - eps_N0 );
      const double sig_N_start = sig_N0 + ( eps_N_start - eps_N0 ) * E_N;
      s_gap_start              = E_N0 * eps_N_start - sig_N_start;
    }
    else if ( s_gap0 <= 0.0 ) { // no transition / transition already finished
      return std::make_pair( sig_N0 + E_N0 * deps_N, E_N0 );
    }
    else { // transition already started
      eps_N_start       = eps_N0;
      deps_N_transition = deps_N;
      s_gap_start       = s_gap0;
    }

    // calculate the strain increment where the transition ends
    const double deps_N_end = -( std::pow( s_gap_start, ( 1.0 - transitionExponent ) ) *
                                 std::pow( transitionWidth, transitionExponent ) ) /
                              ( ( E_N0 - E_N ) * ( 1.0 - transitionExponent ) );

    if ( deps_N_end > deps_N_transition ) { // transition phase ends within this increment
      sig_N = E_N0 * ( eps_N0 + deps_N );
      return std::make_pair( sig_N, E_N0 );
    }

    // integration constant (part of it)
    const double C = std::pow( s_gap_start, ( 1. - transitionExponent ) ) / ( 1. - transitionExponent );

    // s_gap at the end of the step
    const double s_gap = std::pow( ( ( 1.0 - transitionExponent ) *
                                     ( ( E_N0 - E_N ) / std::pow( transitionWidth, transitionExponent ) *
                                         deps_N_transition +
                                       C ) ),
                                   ( 1. / ( 1. - transitionExponent ) ) );

    sig_N = E_N0 * ( eps_N0 + deps_N ) - s_gap;
    return std::make_pair( sig_N, unloadingTransitionStiffness( s_gap, E_N, E_N0 ) );
  }

  double MicroplaneM7Formulation::unloadingTransitionStiffness( const double s_gap,
                                                                const double E_N,
                                                                const double E_N0 ) const
  {
    const double u      = std::min( std::max( s_gap / p.transitionWidth, 0.0 ), 1.0 );
    const double weight = std::pow( u, p.transitionExponent );
    return weight * E_N + ( 1.0 - weight ) * E_N0;
  }

  MicroplaneM7Formulation::NormalResponse MicroplaneM7Formulation::volumetric_boundary(
    const double    sig_v0,
    const Vector6d& strain_macro ) const
  {
    const double E = p.E, k1 = p.k1, k3 = p.k3, k4 = p.k4, k5 = p.k5;

    double E_N0  = normalStiffness(); // microplane normal undamaged stiffness
    double eps_v = ( strain_macro( 0 ) + strain_macro( 1 ) + strain_macro( 2 ) ) / 3;

    Vector6d strain_macroTensorComponents = strain_macro.array() * PInv.array();
    const auto [epsPrinc, depsPrincddeps] = Invariants::principalValuesAndDerivatives( strain_macroTensorComponents );
    Eigen::Index maxIndex, minIndex;
    double       prStrainDiff = epsPrinc.maxCoeff( &maxIndex ) - epsPrinc.minCoeff( &minIndex );

    double alpha = k5 / ( 1 + std::min( std::max( -sig_v0, 0.0 ), c21 ) / E_N0 ) * std::pow( prStrainDiff / k1, c20 ) +
                   k4;
    double sig_vb = -E * k1 * k3 * Math::exp( -eps_v / ( k1 * alpha ) ); // same value for all microplanes

    Vector6d
      princStrainDiffDeriv = ( depsPrincddeps.row( maxIndex ) - depsPrincddeps.row( minIndex ) ).transpose().array() *
                             PInv.array();
    Vector6d tangent = sig_vb * ( -I / ( 3.0 * k1 * alpha ) +
                                  eps_v / ( k1 * alpha * alpha ) * k5 /
                                    ( 1 + std::min( std::max( -sig_v0, 0.0 ), c21 ) / E_N0 ) * c20 *
                                    std::pow( prStrainDiff / k1, c20 - 1. ) / k1 * princStrainDiffDeriv );

    return { sig_vb, tangent };
  }

  MicroplaneM7Formulation::NormalResponse MicroplaneM7Formulation::deviatoric_boundary( const double eps_d,
                                                                                        const double eps_v,
                                                                                        const size_t qp_i ) const
  {
    const double E = p.E, k1 = p.k1;

    double gamma_0 = fc0_reference / E0_reference - p.fc / E;
    double gamma_1 = Math::exp( gamma_0 ) * std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 );
    double beta_2  = c5 * gamma_1 + c7;
    double beta_3  = c6 * gamma_1 + c8;

    double pow2_term = std::max( -eps_d, 0.0 ) / ( k1 * beta_2 );
    double stress    = -E * k1 * beta_3 / ( 1.0 + pow2_term * pow2_term );

    Vector6d tangent = -E * k1 * c6 * Math::exp( gamma_0 ) *
                         ( 1.0 - std::pow( std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 ), 2.0 ) ) * c9 / k1 *
                         Math::heaviside( -eps_v ) * -1.0 / 3.0 * I * 1.0 /
                         ( 1.0 + std::pow( std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ), 2.0 ) ) +
                       E * k1 * beta_3 /
                         std::pow( 1.0 + std::pow( std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ), 2.0 ), 2.0 ) * 2.0 *
                         std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ) *
                         ( Math::heaviside( -eps_d ) / ( k1 * beta_2 ) *
                             ( -ntSplit().pTensors.N.col( qp_i ) + I / 3.0 ) -
                           std::max( -eps_d, 0.0 ) / ( k1 * std::pow( beta_2, 2.0 ) ) * c5 * Math::exp( gamma_0 ) *
                             ( 1.0 - std::pow( std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 ), 2.0 ) ) * c9 / k1 *
                             Math::heaviside( -eps_v ) * -I / 3.0 );

    return { stress, tangent };
  }

  MicroplaneM7Formulation::NormalResponse MicroplaneM7Formulation::tensile_boundary( const double sig_v0,
                                                                                     const double eps_N,
                                                                                     const size_t qp_i ) const
  {
    const double E = p.E, k1 = p.k1, c3 = p.c3;

    double E_N0   = normalStiffness();
    double beta_1 = -c1 + c17 * Math::exp( -c19 * std::max( -sig_v0 - c18, 0.0 ) / E_N0 );

    double stress = E * k1 * beta_1 *
                    Math::exp( -std::max( eps_N - beta_1 * c2 * k1, 0.0 ) /
                               ( c4 * std::max( -sig_v0 / E_N0, 0.0 ) + k1 * c3 ) );
    double tangent = stress * -Math::heaviside( eps_N - beta_1 * c2 * k1 ) /
                     ( c4 * std::max( -sig_v0 / E_N0, 0.0 ) + k1 * c3 );

    return { stress, tangent * ntSplit().pTensors.N.col( qp_i ) };
  }

  MicroplaneM7Formulation::ShearBoundaryInputs MicroplaneM7Formulation::computeShearBoundaryInputs(
    const double    eps_v,
    const double    sig_N,
    const double    sig_M0,
    const double    sig_L0,
    const double    deps_M,
    const double    deps_L,
    const Vector6d& sig_N_tangent_elastic,
    const Vector6d& sig_N_tangent_boundary ) const
  {
    const double k1 = p.k1, k2 = p.k2;

    ShearBoundaryInputs in;
    in.E_T      = shearStiffness(); // there is no damage for shear stiffness
    in.sig_dach = in.E_T * std::max( k1 * c11 - c12 * std::max( eps_v, 0.0 ), 0.0 );

    in.elasticL_prediction  = shear_elasticLaw( sig_L0, deps_L );
    in.elasticM_prediction  = shear_elasticLaw( sig_M0, deps_M );
    in.shear_stress_elastic = std::sqrt( in.elasticL_prediction.stress * in.elasticL_prediction.stress +
                                         in.elasticM_prediction.stress * in.elasticM_prediction.stress );

    in.R27_denom = in.E_T * k1 * k2 + c10 * std::max( in.sig_dach - sig_N, 0.0 );
    in.R27       = c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 / in.R27_denom;

    in.R28_denom = 1. + c10 / ( in.E_T * k1 * k2 ) * in.sig_dach;
    in.R28       = c10 * in.sig_dach / in.R28_denom;

    in.dsigDachddStrain_boundary = in.E_T * Math::heaviside( k1 * c11 - c12 * std::max( eps_v, 0.0 ) ) *
                                   ( -c12 * Math::heaviside( eps_v ) ) * I / 3.0;

    in.dR27ddStrain_elastic = ( c10 * Math::heaviside( in.sig_dach - sig_N ) * ( -sig_N_tangent_elastic ) * in.E_T *
                                k1 * k2 ) /
                                in.R27_denom -
                              c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 /
                                ( in.R27_denom * in.R27_denom ) * c10 * Math::heaviside( in.sig_dach - sig_N ) *
                                ( -sig_N_tangent_elastic );

    in.dR27ddStrain_boundary = ( c10 * Math::heaviside( in.sig_dach - sig_N ) *
                                 ( in.dsigDachddStrain_boundary - sig_N_tangent_boundary ) * in.E_T * k1 * k2 ) /
                                 in.R27_denom -
                               c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 /
                                 ( in.R27_denom * in.R27_denom ) * c10 * Math::heaviside( in.sig_dach - sig_N ) *
                                 ( in.dsigDachddStrain_boundary - sig_N_tangent_boundary );

    in.dR28ddStrain_boundary = c10 * in.dsigDachddStrain_boundary / in.R28_denom -
                               c10 * in.sig_dach / ( in.R28_denom * in.R28_denom ) * ( c10 / ( in.E_T * k1 * k2 ) ) *
                                 in.dsigDachddStrain_boundary;

    return in;
  }

  MicroplaneM7Formulation::ShearResponse MicroplaneM7Formulation::computeShearResponseFromBoundary(
    const ShearBoundaryInputs& in,
    const double               shear_boundary,
    const Vector6d&            shearBoundaryTangent_elastic,
    const Vector6d&            shearBoundaryTangent_boundary,
    const size_t               qp_i ) const
  {
    ShearResponse shear_result;

    if ( in.shear_stress_elastic > shear_boundary && in.shear_stress_elastic != 0.0 ) {
      // boundary value exceeded condition
      shear_result.boundaryExceeded = true;

      shear_result.stressL = in.elasticL_prediction.stress * shear_boundary / in.shear_stress_elastic;
      shear_result.stressM = in.elasticM_prediction.stress * shear_boundary / in.shear_stress_elastic;

      shear_result
        .tangentL_elastic = in.E_T * ntSplit().pTensors.L.col( qp_i ) * shear_boundary / in.shear_stress_elastic +
                            in.elasticL_prediction.stress *
                              ( shearBoundaryTangent_elastic / in.shear_stress_elastic -
                                shear_boundary / ( in.shear_stress_elastic * in.shear_stress_elastic ) * 1. /
                                  in.shear_stress_elastic *
                                  ( in.elasticL_prediction.stress * in.E_T * ntSplit().pTensors.L.col( qp_i ) +
                                    in.elasticM_prediction.stress * in.E_T * ntSplit().pTensors.M.col( qp_i ) ) );

      shear_result.tangentL_boundary = in.elasticL_prediction.stress * shearBoundaryTangent_boundary /
                                       in.shear_stress_elastic;

      shear_result
        .tangentM_elastic = in.E_T * ntSplit().pTensors.M.col( qp_i ) * shear_boundary / in.shear_stress_elastic +
                            in.elasticM_prediction.stress *
                              ( shearBoundaryTangent_elastic / in.shear_stress_elastic -
                                shear_boundary / ( in.shear_stress_elastic * in.shear_stress_elastic ) * 1. /
                                  in.shear_stress_elastic *
                                  ( in.elasticL_prediction.stress * in.E_T * ntSplit().pTensors.L.col( qp_i ) +
                                    in.elasticM_prediction.stress * in.E_T * ntSplit().pTensors.M.col( qp_i ) ) );

      shear_result.tangentM_boundary = in.elasticM_prediction.stress * shearBoundaryTangent_boundary /
                                       in.shear_stress_elastic;
    }
    else {
      // boundary not exceeded -> return elastic values
      shear_result.boundaryExceeded  = false;
      shear_result.stressL           = in.elasticL_prediction.stress;
      shear_result.stressM           = in.elasticM_prediction.stress;
      shear_result.tangentL_elastic  = in.elasticL_prediction.tangent * ntSplit().pTensors.L.col( qp_i );
      shear_result.tangentL_boundary = Vector6d::Zero();
      shear_result.tangentM_elastic  = in.elasticM_prediction.tangent * ntSplit().pTensors.M.col( qp_i );
      shear_result.tangentM_boundary = Vector6d::Zero();
    }

    return shear_result;
  }

  MicroplaneM7Formulation::ShearResponse MicroplaneM7Formulation::shear_law( const double    eps_v,
                                                                             const double    sig_N,
                                                                             const double    sig_M0,
                                                                             const double    sig_L0,
                                                                             const double    deps_M,
                                                                             const double    deps_L,
                                                                             const Vector6d& sig_N_tangent_elastic,
                                                                             const Vector6d& sig_N_tangent_boundary,
                                                                             const size_t    qp_i ) const
  {
    const ShearBoundaryInputs in = computeShearBoundaryInputs( eps_v,
                                                               sig_N,
                                                               sig_M0,
                                                               sig_L0,
                                                               deps_M,
                                                               deps_L,
                                                               sig_N_tangent_elastic,
                                                               sig_N_tangent_boundary );

    const double shear_boundary = sig_N > 0. ? in.R28 : in.R27;

    const Vector6d shearBoundaryTangent_elastic  = sig_N > 0. ? Vector6d::Zero() : in.dR27ddStrain_elastic;
    const Vector6d shearBoundaryTangent_boundary = sig_N > 0. ? in.dR28ddStrain_boundary : in.dR27ddStrain_boundary;

    return computeShearResponseFromBoundary( in,
                                             shear_boundary,
                                             shearBoundaryTangent_elastic,
                                             shearBoundaryTangent_boundary,
                                             qp_i );
  }

  MicroplaneM7Formulation::ElasticPrediction MicroplaneM7Formulation::shear_elasticLaw( const double sig_S0,
                                                                                        const double deps_S ) const
  {
    const double E_T = shearStiffness(); // microplane tangential undamaged stiffness
    return { sig_S0 + E_T * deps_S, E_T };
  }
} // namespace Marmot::Materials
