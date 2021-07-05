#include "Marmot/DryingCreepComplianceSpectrumApproximation.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward.hpp"
#include <cmath>
#include <iostream>
#include <functional>

namespace Marmot::Materials {

    using namespace Eigen;
    using namespace Marmot;

  namespace SpectrumApproximation {

    namespace DryingCreepB3 {


      void computeApproximation( Ref< Solidification::KelvinProperties > kelvinElasticModuli,
                                 Ref< Solidification::KelvinProperties > kelvinRetardationTimes,
                                 const double                            hEnv,
                                 const double                            xi,
                                 const double                            xiZero,
                                 const enum ApproximationOrder           approximationOrder )
      {
        
        const double b = 8. * ( 1. - hEnv );


        switch ( approximationOrder ) {

        case ApproximationOrder::firstOrder:

          auto L1 = [&]( dualDouble tau ) {
            double val = tau.val * autodiff::forward::derivative( phi, autodiff::forward::wrt( tau ), autodiff::forward::at( tau, b, xiZero ) );
            std::cout << val << std::endl;
            return val;
          };

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) *  L1( tau ) );
          }
          break;
        }
      return;
      }

      dualDouble phi( dualDouble xi, double b, double xiZero ){

        return sqrt( f( xi - xiZero, b ) - f(  -xiZero, b  ) );
      }
      
      dualDouble T( dualDouble eta ) { return tanh( eta ); }

      dualDouble S( dualDouble xi ) { return T( dualDouble( sqrt( xi ) ) ); }

      dualDouble f( dualDouble xi, double b ) { return exp( b * S( xi ) ); }

      VectorXd generateRetardationTimes( double min, int num )
      {
        VectorXd kelvinRetardationTimes( num );
        for ( int i = 0; i < num; i++ )
          kelvinRetardationTimes( i ) = min * std::pow( std::sqrt( 10. ), i );
        return kelvinRetardationTimes;
      }
    } // namespace DryingCreepB3
  }   // namespace SpectrumApproximation
} // namespace Marmot::Materials
