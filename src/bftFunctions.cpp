#include "bftFunctions.h"
#include <cmath>

namespace bft{
    namespace Functions
    {
        // return linear interpolation of polynom y at given coordinates (x0, y0) and (x1, y1) at point x 
        double linearInterpolation(double x, double x0, double x1, double y0, double y1)
        { return y0 + (x - x0) * (y1 -y0)/(x1 - x0);}

        // bounded version of std::exp
        double exp(double x)
        {
            if(x <= -64)    // underflow if arg < -708.4 (type double)
                return 0.0;
            if(x >= 64)    // overflow if arg > 709.8 (type double), leave ample margin (e.g. for squaring)
                return std::exp(64);
            return std::exp(x);
        }

    }
}
