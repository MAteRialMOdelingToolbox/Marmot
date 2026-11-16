#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Discretizations::Elements::Spatial1D {
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

    dNdXiSized d2NdXi2( double )
    {
      return dNdXiSized::Zero();
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

    dNdXiSized d2NdXi2( double )
    {
      dNdXiSized d2NdXi2;
      d2NdXi2 << 1., 1., -2.;
      return d2NdXi2;
    }

  } // end of namespace Bar3
} // namespace Marmot::Discretizations::Elements::Spatial1D
