#include "bftMechanics.h"

namespace bft{
    namespace mechanics
    {
        Matrix6 Cel(double E, double nu)
        {
            Matrix6 Cel;
            Cel <<  (1-nu), nu, nu, 0, 0, 0,
                    nu, (1-nu), nu, 0, 0, 0,
                    nu, nu, (1-nu), 0, 0, 0,
                    0, 0, 0, (1-2*nu)/2, 0, 0,
                    0, 0, 0, 0, (1-2*nu)/2, 0,
                    0, 0, 0, 0, 0, (1-2*nu)/2;
            Cel *= E/((1+nu)*(1-2*nu));
            return Cel;
        }

        double macauly(double scalar)
        {
            return scalar >= 0 ? scalar : 0.0;
        }
        int heaviside(double scalar)
        {
            return scalar >= 0 ? 1 : 0;
        }

    }
}
