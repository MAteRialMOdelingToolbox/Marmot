#include "Marmot/MicroplaneM7.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  double MicroplaneM7::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 12 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 12." ) );
    }
    return this->materialProperties[11];
  }

  // constructor
  MicroplaneM7::MicroplaneM7( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      m7_( MicroplaneM7Formulation::readParameters( materialProperties, nMaterialProperties, 7 ) )
  {
    initializeStateLayout();
  }

  void MicroplaneM7::computeStress( state3D&                state,
                                    Marmot::Matrix6d&       dStressDDStrain,
                                    const Marmot::Vector6d& dStrain,
                                    const timeInfo& /*timeInfo*/ ) const
  {
    std::tie( state.stress, dStressDDStrain ) = stress_update( dStrain, state.stateVars );
  }

  std::pair< Vector6d, Matrix6d > MicroplaneM7::stress_update( const Vector6d& dE_macro, double* stateVars ) const
  {
    // get StateVars
    double&    sig_v0        = stateLayout.getAs< double& >( stateVars, "sig_v0" );
    double&    zeta          = stateLayout.getAs< double& >( stateVars, "zeta" );
    mnmp_sized sig_N0        = stateLayout.getAs< mnmp_sized >( stateVars, "sig_N0" );
    mnmp_sized sig_M0        = stateLayout.getAs< mnmp_sized >( stateVars, "sig_M0" );
    mnmp_sized sig_L0        = stateLayout.getAs< mnmp_sized >( stateVars, "sig_L0" );
    mVector6d  strain_macro0 = stateLayout.getAs< mVector6d >( stateVars, "strain_macro0" );
    mnmp_sized eps_N0pos     = stateLayout.getAs< mnmp_sized >( stateVars, "eps_N0pos" );
    mnmp_sized eps_N0neg     = stateLayout.getAs< mnmp_sized >( stateVars, "eps_N0neg" );

    // ------------------------------------------------------------------------------------------------
    // Kinematic projection on the microplanes; colums represent the 6 values for a specific microplane
    // ------------------------------------------------------------------------------------------------
    nmp_sized eps_N  = ( strain_macro0 + dE_macro ).transpose() * m7_.ntSplit().pTensors.N;
    nmp_sized deps_N = dE_macro.transpose() * m7_.ntSplit().pTensors.N;
    nmp_sized deps_M = dE_macro.transpose() * m7_.ntSplit().pTensors.M;
    nmp_sized deps_L = dE_macro.transpose() * m7_.ntSplit().pTensors.L;

    double eps_v0 = ( strain_macro0( 0 ) + strain_macro0( 1 ) + strain_macro0( 2 ) ) /
                    3;                                                     // volumetric strain at start of step
    double deps_v = ( dE_macro( 0 ) + dE_macro( 1 ) + dE_macro( 2 ) ) / 3; // volumetric strain increment

    double    eps_v = eps_v0 + deps_v;                                     // volumetric strain at end of step
    nmp_sized eps_d = ( eps_N.array() - eps_v ).matrix(); // deviatoric microplane strain at end of step

    allmp_response normalResult;
    allmp_response shearLResult;
    allmp_response shearMResult;

    int tensileB_flag = 0, compB_flag = 0, shearB_flag = 0;

    const NormalResponse volumetric_law = m7_.volumetric_boundary( sig_v0,
                                                                   strain_macro0 +
                                                                     dE_macro ); // not microplane dependent

    for ( size_t qp_i = 0; qp_i < nqp; qp_i++ ) {
      // --------------------------------------------------------------------------------------------------------
      // Normal microplane response: elastic prediction vs. tensile/compressive boundaries
      // --------------------------------------------------------------------------------------------------------
      const NormalResponse normal_response = computeMicroplaneNormalResponse( qp_i,
                                                                              sig_N0( qp_i ),
                                                                              deps_N( qp_i ),
                                                                              eps_N( qp_i ),
                                                                              sig_v0,
                                                                              zeta,
                                                                              eps_d( qp_i ),
                                                                              eps_v,
                                                                              volumetric_law,
                                                                              eps_N0pos( qp_i ),
                                                                              eps_N0neg( qp_i ),
                                                                              tensileB_flag,
                                                                              compB_flag );
      normalResult.stress( qp_i )          = normal_response.stress;
      normalResult.tangent.col( qp_i )     = normal_response.tangent;

      // --------------------------------------------------------------------------------------------------------
      // Shear microplane response
      // --------------------------------------------------------------------------------------------------------
      // the boundaries are evaluated with the local strain too, so the full tangent is the sum of both parts
      const auto shear_response = m7_.shear_law( eps_v,
                                                 normalResult.stress( qp_i ),
                                                 sig_M0( qp_i ),
                                                 sig_L0( qp_i ),
                                                 deps_M( qp_i ),
                                                 deps_L( qp_i ),
                                                 normalResult.tangent.col( qp_i ),
                                                 Vector6d::Zero(),
                                                 qp_i );

      shearLResult.stress( qp_i )      = shear_response.stressL;
      shearMResult.stress( qp_i )      = shear_response.stressM;
      shearLResult.tangent.col( qp_i ) = shear_response.tangentL_elastic + shear_response.tangentL_boundary;
      shearMResult.tangent.col( qp_i ) = shear_response.tangentM_elastic + shear_response.tangentM_boundary;
      if ( shear_response.boundaryExceeded )
        shearB_flag++;
    }

    Vector6d stress_macro = m7_.ntSplit().homogenizeStress( normalResult.stress,
                                                            shearMResult.stress,
                                                            shearLResult.stress );

    Matrix6d tangent_macro = m7_.ntSplit().homogenizeTangent( normalResult.tangent,
                                                              shearMResult.tangent,
                                                              shearLResult.tangent );

    // a NaN or inf on any microplane carries through the homogenization
    if ( !stress_macro.allFinite() || !tangent_macro.allFinite() ) {
      MarmotJournal::warningToMSG( MakeString() << "stress:\n"
                                                << stress_macro.transpose() << "\ntangent:\n"
                                                << tangent_macro );
      throw Marmot::StressUpdateFailed( MakeString()
                                        << __PRETTY_FUNCTION__ << ": non-finite value in the stress or the tangent" );
    }

    // commit StateVars
    double sig_v = ( stress_macro( 0 ) + stress_macro( 1 ) + stress_macro( 2 ) ) / 3.;

    zeta += std::max( deps_v, 0.0 );

    sig_v0 = sig_v;
    sig_N0 = normalResult.stress;
    sig_L0 = shearLResult.stress;
    sig_M0 = shearMResult.stress;
    strain_macro0 += dE_macro;
    stateLayout.getAs< double& >( stateVars, "activeTensileBoundary" ) = static_cast< double >( tensileB_flag ) /
                                                                         static_cast< double >( nqp );
    stateLayout.getAs< double& >( stateVars, "activeCompressiveBoundary" ) = static_cast< double >( compB_flag ) /
                                                                             static_cast< double >( nqp );
    stateLayout.getAs< double& >( stateVars, "activeShearBoundary" ) = static_cast< double >( shearB_flag ) /
                                                                       static_cast< double >( nqp );

    return std::make_pair( stress_macro, tangent_macro );
  }

  MicroplaneM7::NormalResponse MicroplaneM7::computeMicroplaneNormalResponse(
    const size_t                        qp_i,
    const double                        sig_N0,
    const double                        deps_N,
    const double                        eps_N,
    const double                        sig_v0,
    const double                        zeta,
    const double                        eps_d,
    const double                        eps_v,
    const MicroplaneM7::NormalResponse& volumetric_law,
    double&                             eps_N0pos,
    double&                             eps_N0neg,
    int&                                tensileB_flag,
    int&                                compB_flag ) const
  {
    // eps_N0pos/eps_N0neg are read here at their pre-increment values (before being
    // conditionally updated below), regardless of which branch ends up governing.
    const NormalResponse
      elasticN_prediction = m7_.normal_elasticLaw( sig_N0, deps_N, eps_N, eps_N0pos, eps_N0neg, sig_v0, zeta, qp_i );

    const NormalResponse tensile_law = m7_.tensile_boundary( sig_v0, eps_N, qp_i );

    if ( elasticN_prediction.stress > tensile_law.stress ) {
      // tensile boundary exceeded
      if ( eps_N > eps_N0pos ) // eps_N0pos shouldn't be negative
        eps_N0pos = eps_N;
      tensileB_flag++;
      return tensile_law;
    }

    const NormalResponse deviatoric_law = m7_.deviatoric_boundary( eps_d, eps_v, qp_i );
    NormalResponse       compressive_boundary;
    compressive_boundary.stress  = deviatoric_law.stress + volumetric_law.stress;
    compressive_boundary.tangent = deviatoric_law.tangent + volumetric_law.tangent;

    if ( elasticN_prediction.stress < compressive_boundary.stress ) {
      // compressive boundary exceeded
      if ( eps_N < eps_N0neg ) // eps_N0neg shouldn't be positive
        eps_N0neg = eps_N;
      compB_flag++;
      return compressive_boundary;
    }
    // elastic behaviour
    return elasticN_prediction;
  }
} // namespace Marmot::Materials
