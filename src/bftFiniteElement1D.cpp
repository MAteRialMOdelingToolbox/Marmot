#include "bftFiniteElement.h"
#include "bftTypedefs.h"

namespace bft {
    namespace FiniteElement {
        namespace Spatial1D {
            namespace Truss2 {
                /*  (0)            (1)
                 *    o------------o
                 *
                 *    <------|-----> xi
                 *     -1    0     1
                 */

                NSized N ( double  xi )
                {
                    NSized N;
                    N <<    ( 1-xi ) * 0.5,
                    ( 1+xi ) * 0.5;
                    return N;
                }

                dNdXiSized dNdXi ( double xi )
                {
                    dNdXiSized dNdXi;
                    dNdXi << -0.5, 0.5;
                    return dNdXi;
                }

            } // end of namespace Truss2


            namespace Truss3 {
                /*
                 * (0)    (2)     (1)
                 *  o------o-------o
                 *
                 *   <-----.-----> xi
                 *   -1    0    1
                 *
                 * */

                NSized N ( double  xi )
                {
                    NSized N;
                    N << ( xi - 1 ) * xi / 2,
                    ( xi + 1 ) * xi / 2,
                    ( 1 - xi*xi );
                    return N;
                }

                dNdXiSized dNdXi ( double xi )
                {
                    dNdXiSized dNdXi;
                    dNdXi <<    xi - 0.5,
                          xi + 0.5,
                          -2*xi;

                    return dNdXi;
                }

            } // end of namespace Truss3
        }
    }
}
