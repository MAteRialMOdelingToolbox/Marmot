#include "Marmot/MarmotViscoelasticity.h"

namespace Marmot::ContinuumMechanics::Viscoelasticity {
  namespace DiscreteSpectrum {

    Properties generateLogarithmicTimes( int n, double min, double spacing )
    {
      Properties times( n );
      for ( int i = 0; i < n; ++i )
        times( i ) = min * std::pow( spacing, i );
      return times;
    }

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta )
    {
      const double dT_tau = dT / tau;
      // Respect extreme values according to Jirasek-Bazant.
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

  } // namespace DiscreteSpectrum
} // namespace Marmot::ContinuumMechanics::Viscoelasticity
