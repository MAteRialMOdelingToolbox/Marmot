#include "Marmot/GM7.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <algorithm>
#include <cmath>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::VoigtNotation;

  namespace {
    // Fallback values for optional material properties 9-12 (c3, fc, transitionWidth,
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

    // Number of substeps used by the explicit (superseded) integration of the unloading
    // transition in transitionFunction().
    constexpr size_t transitionFunctionNumSubsteps = 100;
  } // namespace

  // constructor
  GM7::GM7( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGeneralGradientEnhancedHypoElastic<
        6 >::MarmotMaterialGeneralGradientEnhancedHypoElastic( materialProperties,
                                                               nMaterialProperties,
                                                               materialNumber ),
      // clang-format off
          E( materialProperties[0] ),
          nu( materialProperties[1] ),
          k1( materialProperties[2] ),
          k2( materialProperties[3] ),
          k3( materialProperties[4] ),
          k4( materialProperties[5] ),
          k5( materialProperties[6] ),
          c( materialProperties[7] ),
          m( materialProperties[8] ),
          c3( nMaterialProperties > 9 ? materialProperties[9] : defaultPostpeakTensionSlope ),
          fc( nMaterialProperties > 10 ? materialProperties[10] : defaultCompressiveStrength ),
          transitionWidth( nMaterialProperties > 11 ? materialProperties[11] : defaultTransitionWidth ),
          transitionExponent( nMaterialProperties > 12 ? materialProperties[12] : defaultTransitionExponent )
  // clang-format on
  {
    initializeStateLayout();
  }

  double GM7::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 14 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 14." ) );
    }
    return this->materialProperties[13];
  }

  std::vector< double > GM7::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 15 ) {
      throw std::runtime_error( std::string(
        MakeString() << __PRETTY_FUNCTION__ << ": No nonlocal viscosity given! nMaterialProperties < 15." ) );
    }
    return std::vector< double >( 6, this->materialProperties[14] );
  }

  void GM7::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    Map< const Vector6d > dE_macroLocal( inc.dStrain.data() );
    Map< const Vector6d > dE_macroNonLocal( inc.dK.data() );

    mVector6d       S_macro( res.stress.data() );
    Map< Vector6d > E_macroLocal( res.KLocal.data() );

    mMatrix6d dStressddStrainLocal( tan.dStressddStrain.data() );
    mMatrix6d dStressddStrainNonLocal( tan.dStressddK.data() );
    mMatrix6d dStrainLocalddStrainLocal( tan.dKLocalddStrain.data() );

    // over-nonlocal strain increment used for evaluating the tensile and deviatoric boundaries
    const Vector6d dE_macroOverNonLocal = m * dE_macroNonLocal + ( 1. - m ) * dE_macroLocal;

    const auto [stress_macro, dStressddStrainLocal_, dStressddStrainOverNonLocal] = stress_update( dE_macroLocal,
                                                                                                   dE_macroOverNonLocal,
                                                                                                   res.stateVars );

    // update results
    S_macro                             = stress_macro;
    mVector6d strain_macro0Local        = stateLayout.getAs< mVector6d >( res.stateVars, "strain_macro0Local" );
    mVector6d strain_macro0OverNonLocal = stateLayout.getAs< mVector6d >( res.stateVars, "strain_macro0OverNonLocal" );
    strain_macro0Local += dE_macroLocal;
    strain_macro0OverNonLocal += dE_macroOverNonLocal;
    E_macroLocal = strain_macro0Local;
    res.c.setConstant( c );

    // update tangents
    dStressddStrainLocal      = dStressddStrainLocal_ + dStressddStrainOverNonLocal * ( 1. - m );
    dStressddStrainNonLocal   = dStressddStrainOverNonLocal * m;
    dStrainLocalddStrainLocal = Matrix6d::Identity();
  }

  template < typename T > // so the function accepts Vector6d and Eigen::Map<const Vector6d>
  std::tuple< Vector6d, Matrix6d, Matrix6d > GM7::stress_update( const T&        dE_macroLocal,
                                                                 const Vector6d& dE_macroOverNonLocal,
                                                                 double*         stateVars ) const
  {
    // get StateVars
    double&    sig_v0                    = stateLayout.getAs< double& >( stateVars, "sig_v0" );
    double&    zeta                      = stateLayout.getAs< double& >( stateVars, "zeta" );
    mnmp_sized sig_N0                    = stateLayout.getAs< mnmp_sized >( stateVars, "sig_N0" );
    mnmp_sized sig_M0                    = stateLayout.getAs< mnmp_sized >( stateVars, "sig_M0" );
    mnmp_sized sig_L0                    = stateLayout.getAs< mnmp_sized >( stateVars, "sig_L0" );
    mVector6d  strain_macro0Local        = stateLayout.getAs< mVector6d >( stateVars, "strain_macro0Local" );
    mVector6d  strain_macro0OverNonLocal = stateLayout.getAs< mVector6d >( stateVars, "strain_macro0OverNonLocal" );
    mnmp_sized eps_N0pos                 = stateLayout.getAs< mnmp_sized >( stateVars, "eps_N0pos" );
    mnmp_sized eps_N0neg                 = stateLayout.getAs< mnmp_sized >( stateVars, "eps_N0neg" );

    // ------------------------------------------------------------------------------------------------
    // Kinematic projection on the microplanes; colums represent the 6 values for a specific microplane
    // ------------------------------------------------------------------------------------------------

    // local microplane strains
    nmp_sized eps_N_l = ( strain_macro0Local + dE_macroLocal ).transpose() *
                        ntSplit_.pTensors.N; // eps_N local at the end of the step
    nmp_sized deps_N_l = dE_macroLocal.transpose() * ntSplit_.pTensors.N;
    nmp_sized deps_M_l = dE_macroLocal.transpose() * ntSplit_.pTensors.M;
    nmp_sized deps_L_l = dE_macroLocal.transpose() * ntSplit_.pTensors.L;

    // over non-local microplane strains
    nmp_sized eps_N_onl = ( strain_macro0OverNonLocal + dE_macroOverNonLocal ).transpose() *
                          ntSplit_.pTensors.N; // eps_N over non-local at the end of the step

    double eps_v0_onl = ( strain_macro0OverNonLocal( 0 ) + strain_macro0OverNonLocal( 1 ) +
                          strain_macro0OverNonLocal( 2 ) ) /
                        3;                         // volumetric over non-local strain at start of step
    double deps_v_onl = ( dE_macroOverNonLocal( 0 ) + dE_macroOverNonLocal( 1 ) + dE_macroOverNonLocal( 2 ) ) /
                        3;                         // volumetric over non-local strain increment

    double    eps_v_onl = eps_v0_onl + deps_v_onl; // volumetric over non-local strain at end of step
    nmp_sized eps_d_onl = ( eps_N_onl.array() - eps_v_onl )
                            .matrix();             // deviatoric over non-local microplane strain at end of step

    allmp_response normalResult;
    allmp_response shearLResult;
    allmp_response shearMResult;

    int tensileB_flag = 0, compB_flag = 0, shearB_flag = 0;

    const NormalResponse volumetric_law = volumetric_boundary( sig_v0,
                                                               strain_macro0OverNonLocal +
                                                                 dE_macroOverNonLocal ); // not microplane dependent

    for ( size_t qp_i = 0; qp_i < nqp; qp_i++ ) {
      // --------------------------------------------------------------------------------------------------------
      // Normal microplane response: local-strain elastic prediction vs. over-nonlocal-strain tensile/compressive
      // boundaries
      // --------------------------------------------------------------------------------------------------------
      const NormalResponse normal_response = computeMicroplaneNormalResponse( qp_i,
                                                                              sig_N0( qp_i ),
                                                                              deps_N_l( qp_i ),
                                                                              eps_N_l( qp_i ),
                                                                              eps_N_onl( qp_i ),
                                                                              eps_d_onl( qp_i ),
                                                                              eps_v_onl,
                                                                              sig_v0,
                                                                              zeta,
                                                                              volumetric_law,
                                                                              eps_N0pos( qp_i ),
                                                                              eps_N0neg( qp_i ),
                                                                              tensileB_flag,
                                                                              compB_flag );
      normalResult.stress( qp_i )          = normal_response.stress;
      normalResult.tangent_l.col( qp_i )   = normal_response.tangent_l;
      normalResult.tangent_onl.col( qp_i ) = normal_response.tangent_onl;

      // --------------------------------------------------------------------------------------------------------
      // Shear microplane response
      // --------------------------------------------------------------------------------------------------------
      const ShearResponse shear_response = shear_lawCaner( eps_v_onl,
                                                           normalResult.stress( qp_i ),
                                                           sig_M0( qp_i ),
                                                           sig_L0( qp_i ),
                                                           deps_M_l( qp_i ),
                                                           deps_L_l( qp_i ),
                                                           normalResult.tangent_l.col( qp_i ),
                                                           normalResult.tangent_onl.col( qp_i ),
                                                           qp_i,
                                                           shearB_flag );

      shearLResult.stress( qp_i )          = shear_response.stressL;
      shearLResult.tangent_l.col( qp_i )   = shear_response.tangentL_l;
      shearLResult.tangent_onl.col( qp_i ) = shear_response.tangentL_onl;
      shearMResult.stress( qp_i )          = shear_response.stressM;
      shearMResult.tangent_l.col( qp_i )   = shear_response.tangentM_l;
      shearMResult.tangent_onl.col( qp_i ) = shear_response.tangentM_onl;
    }

    if ( normalResult.stress.array().isNaN().any() ) {
      MarmotJournal::warningToMSG( MakeString() << "stress result normal:\n" << normalResult.stress.transpose() );
      throw Marmot::StressUpdateFailed( MakeString()
                                        << __PRETTY_FUNCTION__ << ": NaN value occurred in normal stress result" );
    }
    if ( normalResult.tangent_l.array().isNaN().any() || normalResult.tangent_onl.array().isNaN().any() ) {
      MarmotJournal::warningToMSG( MakeString() << "normal result tangent (local):\n"
                                                << normalResult.tangent_l.transpose() );
      MarmotJournal::warningToMSG( MakeString() << "normal result tangent (over-nonlocal):\n"
                                                << normalResult.tangent_onl.transpose() );
      throw Marmot::StressUpdateFailed( MakeString()
                                        << __PRETTY_FUNCTION__ << ": NaN value occurred in normal result tangent" );
    }

    Vector6d stress_macro = ntSplit_.homogenizeStress( normalResult.stress, shearMResult.stress, shearLResult.stress );

    Matrix6d dStressddStrainLocal = ntSplit_.homogenizeTangent( normalResult.tangent_l,
                                                                shearMResult.tangent_l,
                                                                shearLResult.tangent_l );

    Matrix6d dStressddStrainOverNonLocal = ntSplit_.homogenizeTangent( normalResult.tangent_onl,
                                                                       shearMResult.tangent_onl,
                                                                       shearLResult.tangent_onl );

    // commit StateVars
    double sig_v = ( stress_macro( 0 ) + stress_macro( 1 ) + stress_macro( 2 ) ) / 3.;

    zeta += std::max( deps_v_onl, 0.0 );

    sig_v0                                                             = sig_v;
    sig_N0                                                             = normalResult.stress;
    sig_M0                                                             = shearMResult.stress;
    sig_L0                                                             = shearLResult.stress;
    stateLayout.getAs< double& >( stateVars, "activeTensileBoundary" ) = static_cast< double >( tensileB_flag ) /
                                                                         static_cast< double >( nqp );
    stateLayout.getAs< double& >( stateVars, "activeCompressiveBoundary" ) = static_cast< double >( compB_flag ) /
                                                                             static_cast< double >( nqp );
    stateLayout.getAs< double& >( stateVars, "activeShearBoundary" ) = static_cast< double >( shearB_flag ) /
                                                                       static_cast< double >( nqp );

    return std::make_tuple( stress_macro, dStressddStrainLocal, dStressddStrainOverNonLocal );
  }

  GM7::NormalResponse GM7::computeMicroplaneNormalResponse( const size_t               qp_i,
                                                            const double               sig_N0,
                                                            const double               deps_N_l,
                                                            const double               eps_N_l,
                                                            const double               eps_N_onl,
                                                            const double               eps_d_onl,
                                                            const double               eps_v_onl,
                                                            const double               sig_v0,
                                                            const double               zeta,
                                                            const GM7::NormalResponse& volumetric_law,
                                                            double&                    eps_N0pos,
                                                            double&                    eps_N0neg,
                                                            int&                       tensileB_flag,
                                                            int&                       compB_flag ) const
  {
    // eps_N0pos/eps_N0neg are read here at their pre-increment values (before being
    // conditionally updated below), regardless of which branch ends up governing.
    const NormalResponse
      elasticN_prediction = normal_elasticLaw( sig_N0, deps_N_l, eps_N_l, eps_N0pos, eps_N0neg, sig_v0, zeta, qp_i );
    const NormalResponse deviatoric_law = deviatoric_boundary( eps_d_onl, eps_v_onl, qp_i );

    NormalResponse compressive_boundary;
    compressive_boundary.stress      = deviatoric_law.stress + volumetric_law.stress;
    compressive_boundary.tangent_l   = deviatoric_law.tangent_l + volumetric_law.tangent_l;
    compressive_boundary.tangent_onl = deviatoric_law.tangent_onl + volumetric_law.tangent_onl;

    const NormalResponse tensile_law = tensile_boundary( sig_v0, eps_N_onl, qp_i );

    if ( elasticN_prediction.stress > tensile_law.stress ) {
      // tensile boundary exceeded
      if ( eps_N_onl > eps_N0pos ) // eps_N0pos shouldn't be negative
        eps_N0pos = eps_N_onl;
      tensileB_flag++;
      return tensile_law;
    }
    if ( elasticN_prediction.stress < compressive_boundary.stress ) {
      // compressive boundary exceeded
      if ( eps_N_onl < eps_N0neg ) // eps_N0neg shouldn't be positive
        eps_N0neg = eps_N_onl;
      compB_flag++;
      return compressive_boundary;
    }
    // elastic behaviour
    return elasticN_prediction;
  }

  GM7::NormalResponse GM7::normal_elasticLaw( const double sig_N0,
                                              const double deps_N,
                                              const double eps_N,
                                              const double eps_N0pos,
                                              const double eps_N0neg,
                                              const double sig_v0,
                                              const double zeta,
                                              const size_t qp_i ) const
  {
    const double E_N0 = E / ( 1.0 - 2.0 * nu ); // microplane normal undamaged stiffness
    double       stress;
    double       tangent;
    const double minStiffness = minStiffnessFraction * E_N0;

    if ( sig_N0 >= 0.0 ) {
      double fzeta       = 1.0 / ( 1.0 + 0.1 * zeta * zeta );
      tangent            = std::max( E_N0 * Math::exp( -c13 * eps_N0pos ) * fzeta, minStiffness );
      double sig_N       = sig_N0 + tangent * deps_N;
      double eps_N0      = eps_N - deps_N;
      double s_gap_trial = E_N0 * eps_N - sig_N;
      if ( s_gap_trial < transitionWidth && ( sig_N0 * deps_N ) < 0.0 &&
           abs( E_N0 - tangent ) > tangentTransitionTol ) {
        std::tie( stress, tangent ) = transitionFunctionClosedForm( eps_N0, deps_N, sig_N0, tangent, E_N0 );
        return { stress, tangent * ntSplit_.pTensors.N.col( qp_i ), Vector6d::Zero() };
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

    return { stress, tangent * ntSplit_.pTensors.N.col( qp_i ), Vector6d::Zero() };
  }

  GM7::NormalResponse GM7::normal_elasticLawUndamaged( const double sig_N0,
                                                       const double deps_N,
                                                       const size_t qp_i ) const
  {
    const double tangent = E / ( 1.0 - 2.0 * nu ); // microplane normal undamaged stiffness
    const double stress  = sig_N0 + tangent * deps_N;
    return { stress, tangent * ntSplit_.pTensors.N.col( qp_i ), Vector6d::Zero() };
  }

  std::pair< double, double > GM7::transitionFunction( const double eps_N0,
                                                       const double deps_N,
                                                       const double sig_N0,
                                                       const double E_N,
                                                       const double E_N0 ) const
  {
    const double s_gap0 = E_N0 * eps_N0 - sig_N0;
    double       eps_N_start, deps_N_transition, sig_N_start;
    if ( s_gap0 > transitionWidth ) { // transition not yet started -> calculate starting point
      eps_N_start       = ( transitionWidth - E_N * eps_N0 + sig_N0 ) / ( E_N0 - E_N );
      deps_N_transition = deps_N - ( eps_N_start - eps_N0 );
      sig_N_start       = sig_N0 + ( eps_N_start - eps_N0 ) * E_N;
    }
    else if ( s_gap0 <= 0.0 ) { // no transition if sig_N0>E_N0*eps_N!
      return std::make_pair( sig_N0 + E_N0 * deps_N, E_N0 );
    }
    else { // transition already started
      eps_N_start       = eps_N0;
      deps_N_transition = deps_N;
      sig_N_start       = sig_N0;
    }

    const size_t n_substeps      = transitionFunctionNumSubsteps;
    double       sigN_transition = sig_N_start;
    double       deps_N_i        = 0.;
    double       ddeps_N         = deps_N_transition / n_substeps;
    double       dsNddeN         = 0.0;

    for ( size_t i = 0; i < n_substeps; i++ ) {
      double s_gap_i = E_N0 * ( eps_N_start + deps_N_i ) - sigN_transition;
      if ( s_gap_i < 0.0 )
        return std::make_pair( E_N0 * ( eps_N0 + deps_N ), E_N0 );

      const auto [E_Nreg, dENregddeNi, dENregddsNi] = unloadingTransitionStiffness( s_gap_i, E_N, E_N0 );

      sigN_transition += E_Nreg * ddeps_N;
      deps_N_i += ddeps_N;

      const double dENregddeN = dENregddeNi * static_cast< double >( i ) / n_substeps + dENregddsNi * dsNddeN;
      dsNddeN += dENregddeN * ddeps_N + E_Nreg / n_substeps;
    }
    return std::make_pair( sigN_transition, dsNddeN );
  }

  std::pair< double, double > GM7::transitionFunctionClosedForm( const double eps_N0,
                                                                 const double deps_N,
                                                                 const double sig_N0,
                                                                 const double E_N,
                                                                 const double E_N0 ) const
  {
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

    sig_N                      = E_N0 * ( eps_N0 + deps_N ) - s_gap;
    const auto [E_Nreg, _, __] = unloadingTransitionStiffness( s_gap, E_N, E_N0 );
    return std::make_pair( sig_N, E_Nreg );
  }

  std::tuple< double, double, double > GM7::unloadingTransitionStiffness( const double s_gap,
                                                                          const double E_N,
                                                                          const double E_N0 ) const
  {
    const double u      = std::min( std::max( s_gap / transitionWidth, 0.0 ), 1.0 );
    const double weight = std::pow( u, transitionExponent );
    const double E_Nreg = weight * E_N + ( 1.0 - weight ) * E_N0;

    // compute tangent
    double dweightddsgap;
    if ( s_gap > 0.0 && s_gap <= transitionWidth )
      dweightddsgap = transitionExponent * std::pow( s_gap / transitionWidth, transitionExponent - 1.0 ) /
                      transitionWidth;
    else
      dweightddsgap = 0.0;

    const double dweigthddeN = dweightddsgap * E_N0;
    const double dweightddsN = -dweightddsgap;
    const double dENregddeN  = dweigthddeN * ( -E_N0 + E_N );
    const double dENregddsN  = dweightddsN * ( -E_N0 + E_N );

    return std::make_tuple( E_Nreg, dENregddeN, dENregddsN );
  }

  GM7::NormalResponse GM7::volumetric_boundary( const double sig_v0, const Vector6d& strain_macro ) const
  {
    double E_N0  = E / ( 1.0 - 2.0 * nu ); // microplane normal undamaged stiffness
    double eps_v = ( strain_macro( 0 ) + strain_macro( 1 ) + strain_macro( 2 ) ) / 3;

    Vector6d strain_macroEngineeringNotation = strain_macro.array() * PInv.array();
    const auto [epsPrinc,
                depsPrincddeps] = Invariants::principalValuesAndDerivatives( strain_macroEngineeringNotation );
    Eigen::Index maxIndex, minIndex;
    double       prStrainDiff = epsPrinc.maxCoeff( &maxIndex ) - epsPrinc.minCoeff( &minIndex );

    double alpha = k5 / ( 1 + std::min( std::max( -sig_v0, 0.0 ), c21 ) / E_N0 ) * std::pow( prStrainDiff / k1, c20 ) +
                   k4;
    double sig_vb = -E * k1 * k3 * exp( -eps_v / ( k1 * alpha ) ); // same value for all microplanes

    Vector6d
      princStrainDiffDeriv = ( depsPrincddeps.row( maxIndex ) - depsPrincddeps.row( minIndex ) ).transpose().array() *
                             PInv.array();
    Vector6d tangent_onl = sig_vb * ( -I / ( 3.0 * k1 * alpha ) +
                                      eps_v / ( k1 * alpha * alpha ) * k5 /
                                        ( 1 + std::min( std::max( -sig_v0, 0.0 ), c21 ) / E_N0 ) * c20 *
                                        std::pow( prStrainDiff / k1, c20 - 1. ) / k1 * princStrainDiffDeriv );

    return { sig_vb, Vector6d::Zero(), tangent_onl };
  }

  GM7::NormalResponse GM7::deviatoric_boundary( const double eps_d, const double eps_v, const size_t qp_i ) const
  {
    double gamma_0 = fc0_reference / E0_reference - fc / E;
    double gamma_1 = Math::exp( gamma_0 ) * std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 );
    double beta_2  = c5 * gamma_1 + c7;
    double beta_3  = c6 * gamma_1 + c8;

    double pow2_term = std::max( -eps_d, 0.0 ) / ( k1 * beta_2 );
    double stress    = -E * k1 * beta_3 / ( 1.0 + pow2_term * pow2_term );

    Vector6d tangent_onl = -E * k1 * c6 * Math::exp( gamma_0 ) *
                             ( 1.0 - std::pow( std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 ), 2.0 ) ) * c9 / k1 *
                             Math::heaviside( -eps_v ) * -1.0 / 3.0 * I * 1.0 /
                             ( 1.0 + std::pow( std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ), 2.0 ) ) +
                           E * k1 * beta_3 /
                             std::pow( 1.0 + std::pow( std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ), 2.0 ), 2.0 ) * 2.0 *
                             std::max( -eps_d, 0.0 ) / ( k1 * beta_2 ) *
                             ( Math::heaviside( -eps_d ) / ( k1 * beta_2 ) *
                                 ( -ntSplit_.pTensors.N.col( qp_i ) + I / 3.0 ) -
                               std::max( -eps_d, 0.0 ) / ( k1 * std::pow( beta_2, 2.0 ) ) * c5 * Math::exp( gamma_0 ) *
                                 ( 1.0 - std::pow( std::tanh( c9 * std::max( -eps_v, 0.0 ) / k1 ), 2.0 ) ) * c9 / k1 *
                                 Math::heaviside( -eps_v ) * -I / 3.0 );

    return { stress, Vector6d::Zero(), tangent_onl };
  }

  GM7::NormalResponse GM7::tensile_boundary( const double sig_v0, const double eps_N, const size_t qp_i ) const
  {
    double E_N0   = E / ( 1.0 - 2.0 * nu );
    double beta_1 = -c1 + c17 * Math::exp( -c19 * std::max( -sig_v0 - c18, 0.0 ) / E_N0 );

    double stress = E * k1 * beta_1 *
                    Math::exp( -std::max( eps_N - beta_1 * c2 * k1, 0.0 ) /
                               ( c4 * std::max( -sig_v0 / E_N0, 0.0 ) + k1 * c3 ) );
    double tangent = stress * -Math::heaviside( eps_N - beta_1 * c2 * k1 ) /
                     ( c4 * std::max( -sig_v0 / E_N0, 0.0 ) + k1 * c3 );

    return { stress, Vector6d::Zero(), tangent * ntSplit_.pTensors.N.col( qp_i ) };
  }

  GM7::ShearBoundaryInputs GM7::computeShearBoundaryInputs( const double    eps_v,
                                                            const double    sig_N,
                                                            const double    sig_M0,
                                                            const double    sig_L0,
                                                            const double    deps_M,
                                                            const double    deps_L,
                                                            const Vector6d& sig_N_tangent_l,
                                                            const Vector6d& sig_N_tangent_onl ) const
  {
    ShearBoundaryInputs in;
    in.E_T = E / ( 1.0 - 2.0 * nu ) * ( 1.0 - 4.0 * nu ) /
             ( 1.0 + nu ); // microplane tangential undamaged stiffness (there is no damage for shear stiffness)
    in.sig_dach = in.E_T * std::max( k1 * c11 - c12 * std::max( eps_v, 0.0 ), 0.0 );

    in.elasticL_prediction  = shear_elasticLaw( sig_L0, deps_L );
    in.elasticM_prediction  = shear_elasticLaw( sig_M0, deps_M );
    in.shear_stress_elastic = std::sqrt( in.elasticL_prediction.stress * in.elasticL_prediction.stress +
                                         in.elasticM_prediction.stress * in.elasticM_prediction.stress );

    in.R27_denom = in.E_T * k1 * k2 + c10 * std::max( in.sig_dach - sig_N, 0.0 );
    in.R27       = c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 / in.R27_denom;

    in.R28_denom = 1. + c10 / ( in.E_T * k1 * k2 ) * in.sig_dach;
    in.R28       = c10 * in.sig_dach / in.R28_denom;

    in.dsigDachddStrain_onl = in.E_T * Math::heaviside( k1 * c11 - c12 * std::max( eps_v, 0.0 ) ) *
                              ( -c12 * Math::heaviside( eps_v ) ) * I / 3.0;

    in.dR27ddStrain_l = ( c10 * Math::heaviside( in.sig_dach - sig_N ) * ( -sig_N_tangent_l ) * in.E_T * k1 * k2 ) /
                          in.R27_denom -
                        c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 /
                          ( in.R27_denom * in.R27_denom ) * c10 * Math::heaviside( in.sig_dach - sig_N ) *
                          ( -sig_N_tangent_l );

    in.dR27ddStrain_onl = ( c10 * Math::heaviside( in.sig_dach - sig_N ) *
                            ( in.dsigDachddStrain_onl - sig_N_tangent_onl ) * in.E_T * k1 * k2 ) /
                            in.R27_denom -
                          c10 * std::max( in.sig_dach - sig_N, 0.0 ) * in.E_T * k1 * k2 /
                            ( in.R27_denom * in.R27_denom ) * c10 * Math::heaviside( in.sig_dach - sig_N ) *
                            ( in.dsigDachddStrain_onl - sig_N_tangent_onl );

    in.dR28ddStrain_onl = c10 * in.dsigDachddStrain_onl / in.R28_denom -
                          c10 * in.sig_dach / ( in.R28_denom * in.R28_denom ) * ( c10 / ( in.E_T * k1 * k2 ) ) *
                            in.dsigDachddStrain_onl;

    if ( std::isnan( in.R27 ) || std::isnan( in.R28 ) ) {
      MarmotJournal::warningToMSG( MakeString() << "R27 = " << in.R27 << ", R28 = " << in.R28 );
      throw Marmot::StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": nan in R27 or R28" );
    }

    return in;
  }

  GM7::ShearResponse GM7::computeShearResponseFromBoundary( const ShearBoundaryInputs& in,
                                                            const double               shear_boundary,
                                                            const Vector6d&            boundary_tangent_l,
                                                            const Vector6d&            boundary_tangent_onl,
                                                            const size_t               qp_i,
                                                            int&                       flag ) const
  {
    ShearResponse shear_result;

    if ( in.shear_stress_elastic > shear_boundary && in.shear_stress_elastic != 0.0 ) {
      // boundary value exceeded condition
      flag++;

      shear_result.stressL = in.elasticL_prediction.stress * shear_boundary / in.shear_stress_elastic;
      shear_result.stressM = in.elasticM_prediction.stress * shear_boundary / in.shear_stress_elastic;

      shear_result.tangentL_l = in.E_T * ntSplit_.pTensors.L.col( qp_i ) * shear_boundary / in.shear_stress_elastic +
                                in.elasticL_prediction.stress *
                                  ( boundary_tangent_l / in.shear_stress_elastic -
                                    shear_boundary / ( in.shear_stress_elastic * in.shear_stress_elastic ) * 1. /
                                      in.shear_stress_elastic *
                                      ( in.elasticL_prediction.stress * in.E_T * ntSplit_.pTensors.L.col( qp_i ) +
                                        in.elasticM_prediction.stress * in.E_T * ntSplit_.pTensors.M.col( qp_i ) ) );

      shear_result.tangentL_onl = in.elasticL_prediction.stress * boundary_tangent_onl / in.shear_stress_elastic;

      shear_result.tangentM_l = in.E_T * ntSplit_.pTensors.M.col( qp_i ) * shear_boundary / in.shear_stress_elastic +
                                in.elasticM_prediction.stress *
                                  ( boundary_tangent_l / in.shear_stress_elastic -
                                    shear_boundary / ( in.shear_stress_elastic * in.shear_stress_elastic ) * 1. /
                                      in.shear_stress_elastic *
                                      ( in.elasticL_prediction.stress * in.E_T * ntSplit_.pTensors.L.col( qp_i ) +
                                        in.elasticM_prediction.stress * in.E_T * ntSplit_.pTensors.M.col( qp_i ) ) );

      shear_result.tangentM_onl = in.elasticM_prediction.stress * boundary_tangent_onl / in.shear_stress_elastic;
    }
    else {
      // boundary not exceeded -> return elastic values
      shear_result.stressL      = in.elasticL_prediction.stress;
      shear_result.stressM      = in.elasticM_prediction.stress;
      shear_result.tangentL_l   = in.elasticL_prediction.tangent * ntSplit_.pTensors.L.col( qp_i );
      shear_result.tangentL_onl = Vector6d::Zero();
      shear_result.tangentM_l   = in.elasticM_prediction.tangent * ntSplit_.pTensors.M.col( qp_i );
      shear_result.tangentM_onl = Vector6d::Zero();
    }

    if ( std::isnan( shear_result.stressL ) || std::isnan( shear_result.stressM ) ) {
      MarmotJournal::warningToMSG( MakeString() << "shear stress L:\n"
                                                << shear_result.stressL << "\nshear stress M:\n"
                                                << shear_result.stressM );
      throw Marmot::StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                                     << ": NaN value occurred in shear stress on microplane " << qp_i );
    }
    if ( shear_result.tangentL_l.array().isNaN().any() || shear_result.tangentL_onl.array().isNaN().any() ||
         shear_result.tangentM_l.array().isNaN().any() || shear_result.tangentM_onl.array().isNaN().any() ) {
      MarmotJournal::warningToMSG( MakeString() << "shear tangent L (local):\n"
                                                << shear_result.tangentL_l << "\nshear tangent M (local):\n"
                                                << shear_result.tangentM_l );
      throw Marmot::StressUpdateFailed(
        MakeString() << __PRETTY_FUNCTION__ << ": NaN value occurred in shear tangent on microplane " << qp_i );
    }

    return shear_result;
  }

  GM7::ShearResponse GM7::shear_lawNguyen( const double    eps_v,
                                           const double    sig_N,
                                           const double    sig_M0,
                                           const double    sig_L0,
                                           const double    deps_M,
                                           const double    deps_L,
                                           const Vector6d& sig_N_tangent_l,
                                           const Vector6d& sig_N_tangent_onl,
                                           const size_t    qp_i,
                                           int&            flag ) const
  {
    const ShearBoundaryInputs in = computeShearBoundaryInputs( eps_v,
                                                               sig_N,
                                                               sig_M0,
                                                               sig_L0,
                                                               deps_M,
                                                               deps_L,
                                                               sig_N_tangent_l,
                                                               sig_N_tangent_onl );

    const double
      s = 0.1 * E * k1; // governs width of transition range (shear boundary under tensile normal stress to compression)

    const double shear_boundary = 0.5 * ( 1 - std::tanh( sig_N / s ) ) * in.R27 +
                                  0.5 * ( 1 + std::tanh( sig_N / s ) ) * in.R28;

    const Vector6d boundary_tangent_l = ( -0.5 ) * ( 1 - std::pow( std::tanh( sig_N / s ), 2.0 ) ) * sig_N_tangent_l /
                                          s * in.R27 +
                                        0.5 * ( 1 - std::tanh( sig_N / s ) ) * in.dR27ddStrain_l +
                                        0.5 * ( 1 - std::pow( std::tanh( sig_N / s ), 2.0 ) ) * sig_N_tangent_l / s *
                                          in.R28;

    const Vector6d boundary_tangent_onl = ( -0.5 ) * ( 1 - std::pow( std::tanh( sig_N / s ), 2.0 ) ) *
                                            sig_N_tangent_onl / s * in.R27 +
                                          0.5 * ( 1 - std::tanh( sig_N / s ) ) * in.dR27ddStrain_onl +
                                          0.5 * ( 1 - std::pow( std::tanh( sig_N / s ), 2.0 ) ) * sig_N_tangent_onl /
                                            s * in.R28 +
                                          0.5 * ( 1 + std::tanh( sig_N / s ) ) * in.dR28ddStrain_onl;

    return computeShearResponseFromBoundary( in, shear_boundary, boundary_tangent_l, boundary_tangent_onl, qp_i, flag );
  }

  GM7::ShearResponse GM7::shear_lawCaner( const double    eps_v,
                                          const double    sig_N,
                                          const double    sig_M0,
                                          const double    sig_L0,
                                          const double    deps_M,
                                          const double    deps_L,
                                          const Vector6d& sig_N_tangent_l,
                                          const Vector6d& sig_N_tangent_onl,
                                          const size_t    qp_i,
                                          int&            flag ) const
  {
    const ShearBoundaryInputs in = computeShearBoundaryInputs( eps_v,
                                                               sig_N,
                                                               sig_M0,
                                                               sig_L0,
                                                               deps_M,
                                                               deps_L,
                                                               sig_N_tangent_l,
                                                               sig_N_tangent_onl );

    const double shear_boundary = sig_N > 0. ? in.R28 : in.R27;

    const Vector6d boundary_tangent_l   = sig_N > 0. ? Vector6d::Zero() : in.dR27ddStrain_l;
    const Vector6d boundary_tangent_onl = sig_N > 0. ? in.dR28ddStrain_onl : in.dR27ddStrain_onl;

    return computeShearResponseFromBoundary( in, shear_boundary, boundary_tangent_l, boundary_tangent_onl, qp_i, flag );
  }

  GM7::ElasticPrediction GM7::shear_elasticLaw( const double sig_S0, const double deps_S ) const
  {
    const double E_T = E / ( 1.0 - 2.0 * nu ) * ( 1.0 - 4.0 * nu ) /
                       ( 1.0 + nu ); // microplane tangential undamaged stiffness
    return { sig_S0 + E_T * deps_S, E_T };
  }
} // namespace Marmot::Materials
