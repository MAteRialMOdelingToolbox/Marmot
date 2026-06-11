#include "Marmot/MarmotWiechert.h"
#include "Marmot/MarmotKelvinChain.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace Wiechert {

    Properties initializeElasticModuli( int nMaxwell, double n )
    {
      Properties elasticModuli( nMaxwell );
      elasticModuli.setConstant( n );
      return elasticModuli;
    }

    Properties initializeRelaxationTimes( int nMaxwell, double m )
    {
      Properties relaxationTimes( nMaxwell );
      relaxationTimes.setConstant( m );
      return relaxationTimes;
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

        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffness += lambda * D * factor;
        dStress += ( 1. - beta ) * stateVars.col( i ).eval() * factor;
      }
    }

    void updateStateVarMatrix( double                dT,
                               const Properties&     elasticModuli,
                               const Properties&     relaxationTimes,
                               Ref< StateVarMatrix > stateVars,
                               const Vector6d&       dStrain,
                               const Matrix6d&       unitD_ijkl )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        KelvinChain::computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars.col( i ) = ( lambda * D ) * unitD_ijkl * dStrain + beta * stateVars.col( i );
      }
    }

  } // namespace Wiechert
} // namespace Marmot::Materials
