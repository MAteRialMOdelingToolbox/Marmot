#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
  namespace FiniteElement {
    namespace Spatial1D {
      namespace Bar2 {
        /*  (0)            (1)
         *    o------------o
         *
         *    <------|-----> xi
         *     -1    0     1
         */

        NSized N( double xi )
        {
          NSized N;
          N << ( 1 - xi ) * 0.5, ( 1 + xi ) * 0.5;
          return N;
        }

        dNdXiSized dNdXi( double xi )
        {
          dNdXiSized dNdXi;
          dNdXi << -0.5, 0.5;
          return dNdXi;
        }

      } // end of namespace Bar2

      namespace Bar3 {
        /*
         * (0)    (2)     (1)
         *  o------o-------o
         *
         *   <-----.-----> xi
         *   -1    0    1
         *
         * */

        NSized N( double xi )
        {
          NSized N;
          N << ( xi - 1 ) * xi / 2, ( xi + 1 ) * xi / 2, ( 1 - xi * xi );
          return N;
        }

        dNdXiSized dNdXi( double xi )
        {
          dNdXiSized dNdXi;
          dNdXi << xi - 0.5, xi + 0.5, -2 * xi;

          return dNdXi;
        }

      } // end of namespace Bar3
    }   // namespace Spatial1D
  }     // namespace FiniteElement
} // namespace Marmot
