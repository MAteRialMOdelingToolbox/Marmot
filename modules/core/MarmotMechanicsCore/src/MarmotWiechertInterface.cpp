#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <ostream>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace WiechertInterface {

    void evaluateWiechert( double                                  dT,
                           const Properties&                       elasticModuli,
                           const Properties&                       relaxationTimes,
                           const StateVarMatrix_force_uu&          stateVars_force_uu,
                           const StateVarMatrix_force_us&          stateVars_force_us,
                           const StateVarMatrix_surface_stress_Z&  stateVars_surface_stress_Z,
                           const StateVarMatrix_surface_stress_Y&  stateVars_surface_stress_Y,
                           const StateVarMatrix_surface_stress_us& stateVars_surface_stress_us,
                           double&                                 uniaxialStiffness,
                           Vector3d&                               dforce_uu,
                           Vector3d&                               dforce_us,
                           Vector9d&                               dsurfaceStress_Z,
                           Vector9d&                               dsurfaceStress_Y,
                           Vector9d&                               dsurfaceStress_us,
                           const double                            factor )
    {
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;

        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffness += lambda * D * factor;
        dforce_uu += ( 1. - beta ) * stateVars_force_uu.col( i ).eval() * factor;
        dforce_us += ( 1. - beta ) * stateVars_force_us.col( i ).eval() * factor;
        dsurfaceStress_Z += ( 1. - beta ) * stateVars_surface_stress_Z.col( i ).eval() * factor;
        dsurfaceStress_Y += ( 1. - beta ) * stateVars_surface_stress_Y.col( i ).eval() * factor;
        dsurfaceStress_us += ( 1. - beta ) * stateVars_surface_stress_us.col( i ).eval() * factor;
      }
    }

    void updateStateVarMatrix_force_uu( double                         dT,
                                        const Properties&              elasticModuli,
                                        const Properties&              relaxationTimes,
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
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_force_uu.col( i ) = ( lambda * D ) * unitH_inv_ij * djumpU + beta * stateVars_force_uu.col( i );
      }
    }

    void updateStateVarMatrix_force_us( double                         dT,
                                        const Properties&              elasticModuli,
                                        const Properties&              relaxationTimes,
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
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_force_us.col( i ) = ( lambda * D ) * unitH_inv_nF_ijk * daverage_strain +
                                      beta * stateVars_force_us.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_Z( double                                 dT,
                                                const Properties&                      elasticModuli,
                                                const Properties&                      relaxationTimes,
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
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_surface_stress_Z.col( i ) = ( lambda * D ) * unitZ_ijkl * dsurfaceStrain +
                                              beta * stateVars_surface_stress_Z.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_Y( double                                 dT,
                                                const Properties&                      elasticModuli,
                                                const Properties&                      relaxationTimes,
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
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_surface_stress_Y.col( i ) = ( lambda * D ) * unitYn_H_inv_Fn_ijkl * dsurfaceStrain +
                                              beta * stateVars_surface_stress_Y.col( i );
      }
    }

    void updateStateVarMatrix_surface_stress_us( double                                  dT,
                                                 const Properties&                       elasticModuli,
                                                 const Properties&                       relaxationTimes,
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
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        auto unit_surface_stress_couple      = djumpU.transpose() * unitH_inv_nF_ijk;
        stateVars_surface_stress_us.col( i ) = ( lambda * D ) * unit_surface_stress_couple.transpose() +
                                               beta * stateVars_surface_stress_us.col( i );
      }
    }

  } // namespace WiechertInterface
} // namespace Marmot::Materials
