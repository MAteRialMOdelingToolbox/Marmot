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
          auto L1 = [&]( dualDouble tau ) {
            double val = tau.val * autodiff::forward::derivative( phi< dualDouble >,
                                                                  autodiff::forward::wrt( tau ),
                                                                  autodiff::forward::at( tau, b, xiZero ) );
            // std::cout << val << std::endl;
            return val;
          };

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * L1( tau ) );
          }
          break;
        }
        case ApproximationOrder::secondOrder: {
          typedef autodiff::forward::HigherOrderDual< 2 > dual2nd;
          const int                                       k = 2;

          auto L2 = [&]( dual2nd tau ) {
            dual2nd tau_ = k * tau;
            double  val_  = -pow( - double( tau_.val ), k ) / factorial( k - 1 );
            val_ *=             autodiff::forward::derivative( phi< dual2nd >,
                                                        autodiff::forward::wrt< k >( tau_ ),
                                                        autodiff::forward::at( tau_, b, xiZero ) );
            // std::cout << val << std::endl;
            return val_;
          };

          // units 1...M
          for ( int i = 0; i < kelvinElasticModuli.size(); i++ ) {
            double tau               = kelvinRetardationTimes( i );
            kelvinElasticModuli( i ) = 1. / ( log( sqrt( 10. ) ) * L2( tau ) );
            kelvinRetardationTimes( i ) *= 0.9 + 0.37* exp( - kelvinRetardationTimes( i ) * kelvinRetardationTimes( i ) / 4 );
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

      int factorial( int n )
      {
        int val = 1;
        for ( int a = 1; a <= n; a++ ) {
          val *= a;
        }
        return val;
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
