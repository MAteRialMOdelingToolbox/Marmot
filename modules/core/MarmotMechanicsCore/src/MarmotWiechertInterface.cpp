#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <ostream>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace WiechertInterface {

    // Properties generateRetardationTimes( int n, double min, double spacing )
    //{ //std::cout<<"Inside generateRetardationTimes"<<std::endl;
    //   Properties retardationTimes( n );
    //   for ( int i = 0; i < n; i++ )
    //     retardationTimes( i ) = min * std::pow( spacing, i );
    //   return retardationTimes;
    // }
    Properties initializeElasticModuli( int nMaxwell, double n )
    {
      Properties elasticModuli( nMaxwell );
      elasticModuli << n;
      return elasticModuli;
    }

    Properties initializeRelaxationTimes( int nMaxwell, double m )
    {
      Properties relaxationTimes( nMaxwell );
      relaxationTimes << m;
      return relaxationTimes;
    }

    void evaluateWiechert( double                           dT,
                           Properties                       elasticModuli,
                           Properties                       relaxationTimes,
                           StateVarMatrix_force_uu          stateVars_force_uu,
                           StateVarMatrix_force_us          stateVars_force_us,
                           StateVarMatrix_surface_stress_Z  stateVars_surface_stress_Z,
                           StateVarMatrix_surface_stress_Y  stateVars_surface_stress_Y,
                           StateVarMatrix_surface_stress_us stateVars_surface_stress_us,
                           double&                          uniaxialStiffness,
                           Vector3d&                        dforce_uu,
                           Vector3d&                        dforce_us,
                           Vector9d&                        dsurfaceStress_Z,
                           Vector9d&                        dsurfaceStress_Y,
                           Vector9d&                        dsurfaceStress_us,
                           const double                     factor )
    {
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;

        computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffness += lambda * D * factor;
        dforce_uu += ( 1. - beta ) * stateVars_force_uu.col( i ).eval() * factor;
        dforce_us += ( 1. - beta ) * stateVars_force_us.col( i ).eval() * factor;
        dsurfaceStress_Z += ( 1. - beta ) * stateVars_surface_stress_Z.col( i ).eval() * factor;
        dsurfaceStress_Y += ( 1. - beta ) * stateVars_surface_stress_Y.col( i ).eval() * factor;
        dsurfaceStress_us += ( 1. - beta ) * stateVars_surface_stress_us.col( i ).eval() * factor;
      }
    }

    void updateStateVarMatrix_force_uu( double                         dT,
                                        Properties                     elasticModuli,
                                        Properties                     relaxationTimes,
                                        Ref< StateVarMatrix_force_uu > stateVars_force_uu,
                                        const Vector3d&                djumpU,
                                        const Matrix3d&                unitH_inv_ij )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_force_uu.col( i ) = ( lambda * D ) * unitH_inv_ij * djumpU + beta * stateVars_force_uu.col( i );
      }
    }

    void updateStateVarMatrix_force_us( double                         dT,
                                        Properties                     elasticModuli,
                                        Properties                     relaxationTimes,
                                        Ref< StateVarMatrix_force_us > stateVars_force_us,
                                        const Vector9d&                daverage_strain,
                                        const Matrix< double, 3, 9 >&  unitH_inv_nF_ijk )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_force_us.col( i ) = ( lambda * D ) * unitH_inv_nF_ijk * daverage_strain +
                                      beta * stateVars_force_us.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_Z( double                                 dT,
                                                Properties                             elasticModuli,
                                                Properties                             relaxationTimes,
                                                Ref< StateVarMatrix_surface_stress_Z > stateVars_surface_stress_Z,
                                                const Vector9d&                        dsurfaceStrain,
                                                const Matrix9d&                        unitZ_ijkl )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_surface_stress_Z.col( i ) = ( lambda * D ) * unitZ_ijkl * dsurfaceStrain +
                                              beta * stateVars_surface_stress_Z.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_Y( double                                 dT,
                                                Properties                             elasticModuli,
                                                Properties                             relaxationTimes,
                                                Ref< StateVarMatrix_surface_stress_Y > stateVars_surface_stress_Y,
                                                const Vector9d&                        dsurfaceStrain,
                                                const Matrix9d&                        unitYn_H_inv_Fn_ijkl )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_surface_stress_Y.col( i ) = ( lambda * D ) * unitYn_H_inv_Fn_ijkl * dsurfaceStrain +
                                              beta * stateVars_surface_stress_Y.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_us( double                                  dT,
                                                 Properties                              elasticModuli,
                                                 Properties                              relaxationTimes,
                                                 Ref< StateVarMatrix_surface_stress_us > stateVars_surface_stress_us,
                                                 const Vector3d&                         djumpU,
                                                 const Matrix< double, 3, 9 >&           unitH_inv_nF_ijk )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        auto unit_surface_stress_couple      = djumpU.transpose() * unitH_inv_nF_ijk;
        stateVars_surface_stress_us.col( i ) = ( lambda * D ) * unit_surface_stress_couple.transpose() +
                                               beta * stateVars_surface_stress_us.col( i );
        // stateVars_surface_stress_us.col( i ) = ( lambda * D ) * unitH_inv_nF_ijk.transpose() * djumpU +
        //                                        beta * stateVars_surface_stress_us.col( i );
      }
    }

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta )
    {
      const double dT_tau = dT / tau;
      // respect extreme values according to Jirasek Bazant
      if ( dT_tau >= 30.0 ) {
        beta   = 0.;
        lambda = 1. / dT_tau;
      }
      else if ( dT_tau < 1e-6 ) {
        beta   = 1.0;
        lambda = 1 - 0.5 * dT_tau + 1. / 6 * dT_tau * dT_tau;
      }
      else {
        beta   = std::exp( -dT_tau );
        lambda = ( 1 - beta ) / dT_tau;
      }
    }

  } // namespace WiechertInterface
} // namespace Marmot::Materials
