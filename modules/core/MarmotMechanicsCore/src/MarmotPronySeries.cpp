#include "Marmot/MarmotPronySeries.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>

namespace Marmot::ContinuumMechanics::Viscoelasticity {
  namespace PronySeries {

    using namespace Marmot;

    void evaluatePronySeries( const Properties&               props,
                              Vector6d&                       stress,
                              Matrix6d&                       stiffness,
                              Eigen::Ref< mapStateVarMatrix > stateVars,
                              const Vector6d&                 dStrain,
                              const double                    dT,
                              const bool                      updateStateVars )
    {
      // ultimate stiffness contribution
      stress += props.ultimateStiffnessMatrix * dStrain;
      stiffness = props.ultimateStiffnessMatrix;

      // prony series terms
      if ( dT > 0.0 ) {
        for ( size_t k = 0; k < props.nPronyTerms; k++ ) {
          const auto& tau       = props.pronyRelaxationTimes.block< 6, 6 >( 0, k * 6 );
          const auto& C         = props.pronyStiffnesses.block< 6, 6 >( 0, k * 6 );
          const auto& currState = stateVars.block< 6, 6 >( 0, k * 6 );

          const Matrix6d eta = tau.cwiseProduct( C );

          const Matrix6d exp_dt_tau = tau.unaryExpr( [dT]( const double& x ) {
            if ( x != 0.0 )
              return std::exp( -dT / x );
            else
              return 0.0;
          } );

          // due to strain increment
          stress += ( eta - eta.cwiseProduct( exp_dt_tau ) ) * dStrain / dT;
          stiffness += ( eta - eta.cwiseProduct( exp_dt_tau ) ) / dT;

          // due to history
          stress -= Matrix6d( currState - exp_dt_tau.cwiseProduct( currState ) ).colwise().sum();

          // update state variables only if it is requested
          if ( updateStateVars )
            stateVars.block< 6, 6 >( 0, 6 * k ) = exp_dt_tau.cwiseProduct( currState ) +
                                                  Matrix6d( eta.array().rowwise() * dStrain.transpose().array() / dT ) -
                                                  Matrix6d( eta.array().rowwise() * dStrain.transpose().array() / dT )
                                                    .cwiseProduct( exp_dt_tau );
        }
      }
    }

    void updateStateVars( const Properties&               props,
                          Eigen::Ref< mapStateVarMatrix > stateVars,
                          const Vector6d&                 dStrain,
                          const double                    dT )
    {
      for ( size_t k = 0; k < props.nPronyTerms; k++ ) {
        const auto& tau       = props.pronyRelaxationTimes.block< 6, 6 >( 0, k * 6 );
        const auto& C         = props.pronyStiffnesses.block< 6, 6 >( 0, k * 6 );
        const auto& currState = stateVars.block< 6, 6 >( 0, k * 6 );

        const Matrix6d eta = tau.cwiseProduct( C );

        const Matrix6d exp_dt_tau = tau.unaryExpr( [dT]( const double& x ) {
          if ( x != 0.0 )
            return std::exp( -dT / x );
          else
            return 0.0;
        } );

        stateVars.block< 6, 6 >( 0, 6 * k ) = exp_dt_tau.cwiseProduct( currState ) +
                                              Matrix6d( eta.array().rowwise() * dStrain.transpose().array() / dT ) -
                                              Matrix6d( eta.array().rowwise() * dStrain.transpose().array() / dT )
                                                .cwiseProduct( exp_dt_tau );
      }
    }

  } // namespace PronySeries
} // namespace Marmot::ContinuumMechanics::Viscoelasticity
