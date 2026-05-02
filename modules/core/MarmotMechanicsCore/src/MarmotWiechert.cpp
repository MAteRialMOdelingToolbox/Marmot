#include "Marmot/MarmotWiechert.h"
#include <iostream>
#include <ostream>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace Wiechert {

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

    void evaluateWiechert( double         dT,
                           Properties     elasticModuli,
                           Properties     relaxationTimes,
                           StateVarMatrix stateVars,
                           double&        uniaxialStiffness,
                           Vector6d&      dStress,
                           const double   factor )
    {
      for ( int i = 0; i < relaxationTimes.size(); i++ ) {
        const double& tau = relaxationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;

        computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffness += lambda * D * factor;
        dStress += ( 1. - beta ) * stateVars.col( i ).eval() * factor;
      }
    }

    void updateStateVarMatrix( double                dT,
                               Properties            elasticModuli,
                               Properties            relaxationTimes,
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
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars.col( i ) = ( lambda * D ) * unitD_ijkl * dStrain + beta * stateVars.col( i );
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

  } // namespace Wiechert
} // namespace Marmot::Materials
