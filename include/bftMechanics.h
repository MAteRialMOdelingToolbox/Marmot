#pragma once 
#include "bftTypedefs.h"

#define isNaN(x) (x!=x)

namespace bft{
    namespace mechanics
    {
        template <typename T> 
            int sgn(T val) {
                return (T(0) < val) - (val < T(0));
            }

        Matrix6 Cel(double E, double nu);
        double macauly(double scalar);
        int heaviside(double scalar);
    }
}
