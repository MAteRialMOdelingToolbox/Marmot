#include "Marmot/MarmotWiechert.h"

namespace Marmot::Materials {

  namespace Wiechert {

    Properties generateRelaxationTimes( int n, double min, double spacing )
    {
      return ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::generateLogarithmicTimes( n, min, spacing );
    }

    void evaluateWiechert( const double                 dT,
                           const Properties&            elasticModuli,
                           const Properties&            relaxationTimes,
                           Eigen::Ref< StateVarMatrix > stateVars,
                           double&                      uniaxialStiffness,
                           Marmot::Vector6d&            dStress,
                           const double                 factor )
    {
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;

        ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffness += lambda * D * factor;
        dStress += ( 1. - beta ) * stateVars.col( i ).eval() * factor;
      }
    }

    void updateStateVarMatrix( const double                 dT,
                               const Properties&            elasticModuli,
                               const Properties&            relaxationTimes,
                               Eigen::Ref< StateVarMatrix > stateVars,
                               const Marmot::Vector6d&      dStrain,
                               const Marmot::Matrix6d&      unitD_ijkl )
    {
      if ( dT <= 1e-14 )
        return;

      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;

        ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars.col( i ) = ( lambda * D ) * unitD_ijkl * dStrain + beta * stateVars.col( i );
      }
    }

  } // namespace Wiechert
} // namespace Marmot::Materials
