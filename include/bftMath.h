#pragma once
#include <string>

namespace bft {
    namespace Math {
        double linearInterpolation ( double x, double x0, double x1, double y0, double y1 );
        double exp ( double x );
        int getExponentPowerTen ( const double x );
        double radToDeg ( const double alpha );
        double degToRad ( const double alpha );
        double macauly ( double scalar );
        int heaviside ( double scalar );
    }
}

