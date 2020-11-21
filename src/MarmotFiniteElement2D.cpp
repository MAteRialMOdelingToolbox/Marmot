#include "MarmotFiniteElement.h"
#include "MarmotTypedefs.h"
#include <iostream>

using namespace Eigen;

namespace Marmot {
    namespace FiniteElement {
        namespace Spatial2D {
            namespace Quad4 {

             NSized N( const Vector2d& xi )
                {
                    /*
                     *   Shape functions
                     *    (3) _______(2)
                     *       |       |
                     *       |       |
                     *       |_______|
                     *    (0)        (1)
                     */

                    NSized N_;
                    // clang-format off
                    N_ <<   1./4 * (1.-xi (0)) * (1.-xi (1)),
                            1./4 * (1.+xi (0)) * (1.-xi (1)),
                            1./4 * (1.+xi (0)) * (1.+xi (1)),
                            1./4 * (1.-xi (0)) * (1.+xi (1));
                    // clang-format on
                    return N_;
                }

                dNdXiSized dNdXi( const Vector2d& xi )
                {

                    dNdXiSized result;

                    // clang-format off
                    result <<   /*          (0)                 (1)                 (2)                 (3) */
                           /* ,xi1 */ -1./4* (1-xi (1)), +1./4* (1-xi (1)),  +1./4* (1+xi (1)),    -1./4* (1+xi (1)), 
                           /*, xi2 */ -1./4* (1-xi (0)), -1./4* (1+xi (0)),  +1./4* (1+xi (0)),    +1./4* (1-xi (0)); 
                    // clang-format on
                    return result;
                }

                Vector2i getBoundaryElementIndices( int faceID )
                {
                    /*
                     *               face 3
                     *           (3) _______(2)
                     *              |       |
                     *       face 4 |       | face 2
                     *              |_______|
                     *           (0)        (1)
                     *               face 1
                     * */
                    switch ( faceID ) {
                    case 1: {
                        return ( Vector2i() << 0, 1 ).finished();
                    }
                    case 2: {
                        return ( Vector2i() << 1, 2 ).finished();
                    }
                    case 3: {
                        return ( Vector2i() << 2, 3 ).finished();
                    }
                    case 4: {
                        return ( Vector2i() << 3, 0 ).finished();
                    }
                    default: {
                        throw std::invalid_argument( "Quad4: invalid face ID specifed" );
                    }
                    }
                }


            } // end of namespace Quad4
            namespace Quad8 {
                NSized N( const Vector2d& xi )
                {
                    /* Shape functions
                     *
                     *         (6)
                     *   (3) _______(2)
                     *      |       |
                     *  (7) |       | (5)
                     *      |_______|
                     *   (0)        (1)
                     *         (4)
                     *
                     * */

                    const double              xi0 = xi( 0 );
                    const double              xi1 = xi( 1 );
                    Matrix<double, nNodes, 1> N_;
                    // clang-format off
                    N_ <<   (1.-xi0 ) * (1.-xi1 ) * (-1-xi0-xi1 ) /4,
                            (1.+xi0 ) * (1.-xi1 ) * (-1+xi0-xi1 ) /4,
                            (1.+xi0 ) * (1.+xi1 ) * (-1+xi0+xi1 ) /4,
                            (1.-xi0 ) * (1.+xi1 ) * (-1-xi0+xi1 ) /4,
                            (1-xi0*xi0 )       * (1-xi1      ) /2,
                            (1+xi0      )     * (1-xi1*xi1 ) /2,
                            (1-xi0*xi0 )       * (1+xi1     ) /2,
                            (1-xi0      )     * (1-xi1*xi1 ) /2;

                    // clang-format on
                    return N_;
                }

                dNdXiSized dNdXi( const Vector2d& xi )
                {

                    Matrix<double, nDim, nNodes> result;

                    const double xi0 = xi( 0 );
                    const double xi1 = xi( 1 );
                    // clang-format off
                    result <<
                           /*                  1                           2                           3                           4
                            *                  5                           6                           7                           8 */
                           /* ,Xi1 */     (1-xi1) * (2*xi0+xi1) /4,  (1-xi1) * (2*xi0-xi1) /4, (1+xi1) * (2*xi0+xi1) /4, (1+xi1) * (2*xi0-xi1) /4,
                           /* ,Xi1 */     -xi0* (1-xi1),               (1-xi1*xi1) /2,              -xi0* (1+xi1),               - (1-xi1*xi1) /2,
                           /*                  1                           2                           3                           4
                            *                  5                           6                           7                           8*/
                           /* ,Xi2 */     (1-xi0) * (xi0+2*xi1) /4,      - (1+xi0) * (xi0-2*xi1) /4,     (1+xi0) * (xi0+2*xi1) /4,      (1-xi0) * (-xi0+2*xi1) /4,
                           /* ,Xi2 */     - (1-xi0*xi0) /2,             -xi1* (1+xi0),               + (1-xi0*xi0) /2,             -xi1* (1-xi0);
                    // clang-format on  

                    return result;
                }

                Vector3i getBoundaryElementIndices ( int faceID )
                {
                    switch ( faceID ) {
                    case 1: {
                        return ( Vector3i() << 0,1,4 ).finished();
                    }
                    case 2: {
                        return ( Vector3i() << 1,2,5 ).finished();
                    }
                    case 3: {
                        return ( Vector3i() << 2,3,6 ).finished();
                    }
                    case 4: {
                        return ( Vector3i() << 3,0,7 ).finished();
                    }
                    default: {
                        throw std::invalid_argument ( "Quad8: invalid face ID specifed" );
                    }
                    }
                }

            } // end of namespace Quad8

            void modifyCharElemLengthAbaqusLike ( double& charElemLength, int intPoint )
            {
                switch ( intPoint ) {
                case 0:         // central node
                    charElemLength *= 2./3.;
                    break;
                case 1:         // corner nodes
                case 2:
                case 3:
                case 4:
                    charElemLength *= 5./12;
                    break;
                case 5:         // middle nodes
                case 6:
                case 7:
                case 8:
                    charElemLength *= std::sqrt ( 5./18. );
                    break;
                }
            }
        } // end of namespace Spatial2D
    } // end of namespace FiniteElement
} // end of namespace Marmot
