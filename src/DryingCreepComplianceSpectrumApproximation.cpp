#include "Marmot/DryingCreepComplianceSpectrumApproximation.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward.hpp"
#include <autodiff/forward/forward.hpp>
#include <cmath>
#include <functional>
#include <iostream>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;
  using namespace std;

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

        case ApproximationOrder::firstOrder: {

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * computeSpectrum< 1 >( tau, b, xiZero ) );
          }
          break;
        }
        case ApproximationOrder::secondOrder: {

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * computeSpectrum< 2 >( tau, b, xiZero ) );
            kelvinRetardationTimes( i ) *= 0.9 +
                                           0.37 * exp( -kelvinRetardationTimes( i ) * kelvinRetardationTimes( i ) / 4 );
          }
          break;
        }
        case ApproximationOrder::thirdOrder: {

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * computeSpectrum< 3 >( tau, b, xiZero ) );
            kelvinRetardationTimes( i ) *= 0.9 +
                                           0.37 * exp( -kelvinRetardationTimes( i ) * kelvinRetardationTimes( i ) / 4 );
          }
          break;
        }
          /*
          case ApproximationOrder::fourthOrder: {
            typedef autodiff::forward::HigherOrderDual< 4 > dual4th;
            const int                                       k = 4;

            auto L4 = [&]( dual4th tau ) {
              dual4th tau__ = k * tau;
              double  val__  = -pow( - double( tau__.val ), k ) / factorial( k - 1 );
              val__ *=             autodiff::forward::derivative( phi< dual4th >,
                                                          autodiff::forward::wrt< k >( tau__ ),
                                                          autodiff::forward::at( tau__, b, xiZero ) );
              // std::cout << val << std::endl;
              return val__;
            };

            // units 1...M
            for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
              double tau               = kelvinRetardationTimes( i );
              kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * L4( tau ) );
            }
            break;
          }
          */
        }
        return;
      }

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
