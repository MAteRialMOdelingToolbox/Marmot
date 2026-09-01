#include "Marmot/MarmotKelvinChain.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace KelvinChain {

    Properties generateRetardationTimes( int n, double min, double spacing )
    {
      return ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::generateLogarithmicTimes( n, min, spacing );
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
      ContinuumMechanics::Viscoelasticity::DiscreteSpectrum::computeLambdaAndBeta( dT, tau, lambda, beta );
    }

  } // namespace KelvinChain
} // namespace Marmot::Materials
