#include "MarmotConstants.h"
#include "MarmotTypedefs.h"
#include <cmath>

namespace Marmot {
    namespace Constants {
        double cubicRootEps()
        {
            const static double cubicRootEps = std::pow( std::numeric_limits<double>::epsilon(), 1. / 3 );
            return cubicRootEps;
        }
        double squareRootEps()
        {
            const static double squareRootEps = std::pow( std::numeric_limits<double>::epsilon(), 1. / 2 );
            return squareRootEps;
        }

        const double SquareRootEps = squareRootEps();
        const double CubicRootEps = cubicRootEps();
    } // namespace Constants
} // namespace Marmot
