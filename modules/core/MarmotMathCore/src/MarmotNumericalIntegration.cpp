#include "Marmot/MarmotNumericalIntegration.h"
#include <iostream>
#include <iterator>

namespace Marmot {
  namespace NumericalAlgorithms::Integration {

    double integrateScalarFunction( scalar_to_scalar_function_type     f,
                                    const std::tuple< double, double > integrationLimits,
                                    const int                          n,
                                    const integrationRule              intRule )
    {

      // create linear spacing
      double deltaX = ( std::get< 1 >( integrationLimits ) - std::get< 0 >( integrationLimits ) ) / ( n );
      std::vector< double > xValues( n + 1 );

      std::generate( xValues.begin(), xValues.end(), [n = 0, &deltaX]() mutable { return n++ * deltaX; } );
      double val = 0.;

      switch ( intRule ) {
      case integrationRule::midpoint:
        for ( int i = 0; i < n; i++ ) {
          val += f( ( xValues[i + 1] + xValues[i] ) / 2. ) * deltaX;
        }
        break;
      case integrationRule::trapezodial:
        for ( int i = 0; i < n; i++ ) {
          val += ( f( xValues[i + 1] ) + f( xValues[i] ) ) / 2. * deltaX;
        }
        break;
      case integrationRule::simpson:
        for ( int i = 0; i < n; i++ ) {
          val += deltaX / 6. *
                 ( f( xValues[i] ) + 4 * f( ( xValues[i + 1] + xValues[i] ) / 2. ) + f( xValues[i + 1] ) );
        }
        break;
      default: throw std::invalid_argument( "Invalid integration rule!" );
      }
      return val;
    }
  } // namespace NumericalAlgorithms::Integration
} // namespace Marmot
