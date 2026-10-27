#include "Marmot/GM7.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  // constructor
  GM7::GM7( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGeneralGradientEnhancedHypoElastic<
        6 >::MarmotMaterialGeneralGradientEnhancedHypoElastic( materialProperties,
                                                               nMaterialProperties,
                                                               materialNumber ),
      m7_( MicroplaneM7Formulation::readParameters( materialProperties, nMaterialProperties, 9 ) ),
      c( materialProperties[7] ),
      m( materialProperties[8] )
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
    const Vector6d&       dE_macroLocal = inc.dStrain;
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

  std::tuple< Vector6d, Matrix6d, Matrix6d > GM7::stress_update( const Vector6d& dE_macroLocal,
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
                        m7_.ntSplit.pTensors.N; // eps_N local at the end of the step
    nmp_sized deps_N_l = dE_macroLocal.transpose() * m7_.ntSplit.pTensors.N;
    nmp_sized deps_M_l = dE_macroLocal.transpose() * m7_.ntSplit.pTensors.M;
    nmp_sized deps_L_l = dE_macroLocal.transpose() * m7_.ntSplit.pTensors.L;

    // over non-local microplane strains
    nmp_sized eps_N_onl = ( strain_macro0OverNonLocal + dE_macroOverNonLocal ).transpose() *
                          m7_.ntSplit.pTensors.N; // eps_N over non-local at the end of the step

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

    const auto volumetric_law = m7_.volumetric_boundary( sig_v0,
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
      const auto shear_response = m7_.shear_law( eps_v_onl,
                                                 normalResult.stress( qp_i ),
                                                 sig_M0( qp_i ),
                                                 sig_L0( qp_i ),
                                                 deps_M_l( qp_i ),
                                                 deps_L_l( qp_i ),
                                                 normalResult.tangent_l.col( qp_i ),
                                                 normalResult.tangent_onl.col( qp_i ),
                                                 qp_i );
      if ( shear_response.boundaryExceeded )
        shearB_flag++;

      // the elastic predictions are evaluated with the local, the boundaries with the over-nonlocal strain
      shearLResult.stress( qp_i )          = shear_response.stressL;
      shearLResult.tangent_l.col( qp_i )   = shear_response.tangentL_elastic;
      shearLResult.tangent_onl.col( qp_i ) = shear_response.tangentL_boundary;
      shearMResult.stress( qp_i )          = shear_response.stressM;
      shearMResult.tangent_l.col( qp_i )   = shear_response.tangentM_elastic;
      shearMResult.tangent_onl.col( qp_i ) = shear_response.tangentM_boundary;
    }

    Vector6d stress_macro = m7_.ntSplit.homogenizeStress( normalResult.stress,
                                                          shearMResult.stress,
                                                          shearLResult.stress );

    Matrix6d dStressddStrainLocal = m7_.ntSplit.homogenizeTangent( normalResult.tangent_l,
                                                                   shearMResult.tangent_l,
                                                                   shearLResult.tangent_l );

    Matrix6d dStressddStrainOverNonLocal = m7_.ntSplit.homogenizeTangent( normalResult.tangent_onl,
                                                                          shearMResult.tangent_onl,
                                                                          shearLResult.tangent_onl );

    // a NaN or inf on any microplane carries through the homogenization
    if ( !stress_macro.allFinite() || !dStressddStrainLocal.allFinite() || !dStressddStrainOverNonLocal.allFinite() ) {
      MarmotJournal::warningToMSG( MakeString() << "stress:\n"
                                                << stress_macro.transpose() << "\ntangent (local):\n"
                                                << dStressddStrainLocal << "\ntangent (over-nonlocal):\n"
                                                << dStressddStrainOverNonLocal );
      throw Marmot::StressUpdateFailed( MakeString()
                                        << __PRETTY_FUNCTION__ << ": non-finite value in the stress or the tangents" );
    }

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

  GM7::NormalResponse GM7::computeMicroplaneNormalResponse(
    const size_t                                   qp_i,
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
    int&                                           compB_flag ) const
  {
    // eps_N0pos/eps_N0neg are read here at their pre-increment values (before being
    // conditionally updated below), regardless of which branch ends up governing.
    // The elastic prediction depends on the local strain only, the boundaries on the over-nonlocal strain only.
    const auto elasticN_prediction = m7_.normal_elasticLaw( sig_N0,
                                                            deps_N_l,
                                                            eps_N_l,
                                                            eps_N0pos,
                                                            eps_N0neg,
                                                            sig_v0,
                                                            zeta,
                                                            qp_i );
    const auto tensile_law         = m7_.tensile_boundary( sig_v0, eps_N_onl, qp_i );

    if ( elasticN_prediction.stress > tensile_law.stress ) {
      // tensile boundary exceeded
      if ( eps_N_onl > eps_N0pos ) // eps_N0pos shouldn't be negative
        eps_N0pos = eps_N_onl;
      tensileB_flag++;
      return { tensile_law.stress, Vector6d::Zero(), tensile_law.tangent };
    }

    const auto   deviatoric_law       = m7_.deviatoric_boundary( eps_d_onl, eps_v_onl, qp_i );
    const double compressive_boundary = deviatoric_law.stress + volumetric_law.stress;
    if ( elasticN_prediction.stress < compressive_boundary ) {
      // compressive boundary exceeded
      if ( eps_N_onl < eps_N0neg ) // eps_N0neg shouldn't be positive
        eps_N0neg = eps_N_onl;
      compB_flag++;
      return { compressive_boundary, Vector6d::Zero(), deviatoric_law.tangent + volumetric_law.tangent };
    }
    // elastic behaviour
    return { elasticN_prediction.stress, elasticN_prediction.tangent, Vector6d::Zero() };
  }
} // namespace Marmot::Materials
