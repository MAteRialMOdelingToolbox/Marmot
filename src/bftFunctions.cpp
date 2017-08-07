#include "bftFunctions.h"
#include <cmath>

namespace bft{
    namespace Functions
    {
        // return linear interpolation of polynom y at given coordinates (x0, y0) and (x1, y1) at point x 
        double linearInterpolation(double x, double x0, double x1, double y0, double y1)
        { return y0 + (x - x0) * (y1 -y0)/(x1 - x0);}   

        double exp(double x)
        {
            // bounded version of std::exp
            if(x <= -708.4)
                return 0.0;
            if(x >= 709.8)
                return std::exp(709.8);
            return std::exp(x);
        }

    }
}
