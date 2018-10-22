#include "bftConstants.h"
#include "bftMath.h"
#include <cmath>
#include <math.h>

namespace bft {
    namespace Math {
        // return linear interpolation of polynom y at given coordinates (x0, y0) and (x1, y1) at
        // point x
        double linearInterpolation( double x, double x0, double x1, double y0, double y1 ) {
            return y0 + ( x - x0 ) * ( y1 - y0 ) / ( x1 - x0 );
        }

        // bounded version of std::exp
        double exp( double x ) {
            if ( x <= -64 ) // underflow if arg < -708.4 (type double)
                return 0.0;
            if ( x >= 64 ) // overflow if arg > 709.8 (type double), leave ample margin (e.g. for
                           // squaring)
                return std::exp( 64 );
            return std::exp( x );
        }

        // return the exponent to the power of ten of an expression like 5*10^5 --> return 5
        int getExponentPowerTen( const double x ) {
            if ( x >= 1e-16 ) // positive number
                return floor( log10( x ) );
            else if ( x <= 1e-16 ) // negative number
                return floor( log10( abs( x ) ) );
            else // number close to 0
                return 0;
        }

        double radToDeg( const double alpha ) { return alpha * 180 / bft::Constants::Pi; }

        double degToRad( const double alpha ) { return alpha / 180 * bft::Constants::Pi; }

        double macauly( double scalar ) { return scalar >= 0 ? scalar : 0.0; }

        int heaviside( double scalar ) { return scalar >= 0 ? 1 : 0; }

    } // namespace Math
} // namespace bft
