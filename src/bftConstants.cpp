#include "bftConstants.h"
#include "bftTypedefs.h"
#include <cmath>

namespace bft{
    namespace Constants
    {
        double cubicRootEps()
            {const static double cubicRootEps = std::pow(std::numeric_limits<double>::epsilon(), 1./3);
                return cubicRootEps;}
        double squareRootEps()
            {const static double squareRootEps = std::pow(std::numeric_limits<double>::epsilon(), 1./2);
                return squareRootEps;}
    }
    namespace Functions
    {
        // return linear interpolation of polynom y at given coordinates (x0, y0) and (x1, y1) at point x 
        double linearInterpolation(double x, double x0, double x1, double y0, double y1)
        { return y0 + (x - x0) * (y1 -y0)/(x1 - x0);}   

    }
}
