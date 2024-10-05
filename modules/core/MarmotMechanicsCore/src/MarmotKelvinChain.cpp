#include "Marmot/MarmotKelvinChain.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace KelvinChain {

    Properties generateRetardationTimes( int n, double min, double spacing )
    {
      Properties retardationTimes( n );
      for ( int i = 0; i < n; i++ )
        retardationTimes( i ) = min * std::pow( spacing, i );
      return retardationTimes;
    }

    void evaluateKelvinChain( double         dT,
                              Properties     elasticModuli,
                              Properties     retardationTimes,
                              StateVarMatrix stateVars,
                              double&        uniaxialCompliance,
                              Vector6d&      dStrain,
                              const double   factor )
    {

      for ( int i = 0; i < retardationTimes.size(); i++ ) {
        const double& tau = retardationTimes( i );
        const double& D   = elasticModuli( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialCompliance += ( 1. - lambda ) / D * factor;
        dStrain += ( 1. - beta ) * stateVars.col( i ) * factor;
      }
    }

    void updateStateVarMatrix( double                dT,
                               Properties            elasticModuli,
                               Properties            retardationTimes,
                               Ref< StateVarMatrix > stateVars,
                               const Vector6d&       dStress,
                               const Matrix6d&       unitComplianceMatrix )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < retardationTimes.size(); i++ ) {
        const double& tau = retardationTimes( i );
        const double& D   = elasticModuli( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars.col( i ) = ( lambda / D ) * unitComplianceMatrix * dStress + beta * stateVars.col( i );
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

  } // namespace KelvinChain
} // namespace Marmot::Materials
