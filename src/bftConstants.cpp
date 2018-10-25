#include "bftConstants.h"
#include "bftTypedefs.h"
#include <cmath>

namespace bft {
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
    } // namespace Constants
} // namespace bft
