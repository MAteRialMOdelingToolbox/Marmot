#include "bftFiniteElement.h"
#include "bftTypedefs.h"
#include <iostream>

namespace bft {
    namespace FiniteElement {
        namespace Spatial3D {
            namespace Hexa8 {
                NSized N( const Vector3d& xi ) {
                    /*
                     *        (8)________(7)
                     *          /|      /|
                     *      (5)/_|_____(6)
                     *        |  |     | |
                     *        |(4)_____|_|(3)
                     *        | /      | /
                     *        |/_______|/
                     *     (1)        (2)
                     *
                     *
                     *     x3  x2
                     *     | /
                     *     |/___x1
                     *
                     *  */

                    NSized N_;

                    // clang-format off
                    N_ << (1-xi (0)) * (1-xi (1)) * (1-xi (2)),
                    (1+xi (0)) * (1-xi (1)) * (1-xi (2)),
                    (1+xi (0)) * (1+xi (1)) * (1-xi (2)),
                    (1-xi (0)) * (1+xi (1)) * (1-xi (2)),
                    (1-xi (0)) * (1-xi (1)) * (1+xi (2)),
                    (1+xi (0)) * (1-xi (1)) * (1+xi (2)),
                    (1+xi (0)) * (1+xi (1)) * (1+xi (2)),
                    (1-xi (0)) * (1+xi (1)) * (1+xi (2));
                    // clang-format on  

                    N_ *= 1./8;
                    return N_;

                }
                dNdXiSized dNdXi ( const Vector3d& xi )
                {
                    dNdXiSized dNdXi_;

                    // clang-format off
                    dNdXi_ <<
                           -1 * (1-xi (1)) * (1-xi (2)),
                           +1 * (1-xi (1)) * (1-xi (2)),
                           +1 * (1+xi (1)) * (1-xi (2)),
                           -1 * (1+xi (1)) * (1-xi (2)),
                           -1 * (1-xi (1)) * (1+xi (2)),
                           +1 * (1-xi (1)) * (1+xi (2)),
                           +1 * (1+xi (1)) * (1+xi (2)),
                           -1 * (1+xi (1)) * (1+xi (2)),

                           (1-xi (0)) * -1 * (1-xi (2)),
                           (1+xi (0)) * -1 * (1-xi (2)),
                           (1+xi (0)) * +1 * (1-xi (2)),
                           (1-xi (0)) * +1 * (1-xi (2)),
                           (1-xi (0)) * -1 * (1+xi (2)),
                           (1+xi (0)) * -1 * (1+xi (2)),
                           (1+xi (0)) * +1 * (1+xi (2)),
                           (1-xi (0)) * +1 * (1+xi (2)),

                           (1-xi (0)) * (1-xi (1)) * -1,
                           (1+xi (0)) * (1-xi (1)) * -1,
                           (1+xi (0)) * (1+xi (1)) * -1,
                           (1-xi (0)) * (1+xi (1)) * -1,
                           (1-xi (0)) * (1-xi (1)) * +1,
                           (1+xi (0)) * (1-xi (1)) * +1,
                           (1+xi (0)) * (1+xi (1)) * +1,
                           (1-xi (0)) * (1+xi (1)) * +1;

                    dNdXi_ *= 1./8;
                    // clang-format on

                    return dNdXi_;
                }

                Matrix3d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates ) {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes>( dNdXi, coordinates );
                }

                BSized B( const dNdXiSized& dNdX ) {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial3D::B<nNodes>( dNdX );
                }

                Vector4i getBoundaryElementIndices( int faceID ) {
                    switch ( faceID ) {
                    case 1: {
                        return ( Vector4i() << 3, 2, 1, 0 ).finished();
                    }
                    case 2: {
                        return ( Vector4i() << 4, 5, 6, 7 ).finished();
                    }
                    case 3: {
                        return ( Vector4i() << 0, 1, 5, 4 ).finished();
                    }
                    case 4: {
                        return ( Vector4i() << 6, 5, 1, 2 ).finished();
                    }
                    case 5: {
                        return ( Vector4i() << 7, 6, 2, 3 ).finished();
                    }
                    case 6: {
                        return ( Vector4i() << 4, 7, 3, 0 ).finished();
                    }
                    default: { throw std::invalid_argument( "Hexa8: invalid face ID specifed" ); }
                    }
                }

            } // end of namespace Hexa8

            namespace Hexa20 {
                NSized N( const Vector3d& xi ) {
                    /*
                     *           (8)_____(15)_____(7)
                     *            /|             /|
                     *       (16)/ |            /(14)
                     *          /  |           /  |
                     *      (5)/___|__(13)___(6)  |
                     *        |  (18)         | (17)
                     *        |    |          |   |
                     *        |  (4)___(11)___|___|(3)
                     *      (15)  /         (16)  /
                     *        |  /(12)        |  /(10)
                     *        | /             | /
                     *        |/______________|/
                     *     (1)       (9)      (2)
                     *
                     *
                     *     x3  x2
                     *     | /
                     *     |/___x1
                     *
                     *  */

                    NSized N_;
                    // Taken from mpFEM - Peter Gamnitzer
                    // clang-format off
                    N_ <<
                       -0.125* (1-xi (0)) * (1-xi (1)) * (1-xi (2)) * (2+xi (0)+xi (1)+xi (2)),
                       -0.125* (1+xi (0)) * (1-xi (1)) * (1-xi (2)) * (2-xi (0)+xi (1)+xi (2)),
                       -0.125* (1+xi (0)) * (1+xi (1)) * (1-xi (2)) * (2-xi (0)-xi (1)+xi (2)),
                       -0.125* (1-xi (0)) * (1+xi (1)) * (1-xi (2)) * (2+xi (0)-xi (1)+xi (2)),
                       -0.125* (1-xi (0)) * (1-xi (1)) * (1+xi (2)) * (2+xi (0)+xi (1)-xi (2)),
                       -0.125* (1+xi (0)) * (1-xi (1)) * (1+xi (2)) * (2-xi (0)+xi (1)-xi (2)),
                       -0.125* (1+xi (0)) * (1+xi (1)) * (1+xi (2)) * (2-xi (0)-xi (1)-xi (2)),
                       -0.125* (1-xi (0)) * (1+xi (1)) * (1+xi (2)) * (2+xi (0)-xi (1)-xi (2)),
                       0.25* (1-xi (0) *xi (0)) * (1-xi (1)) * (1-xi (2)),
                       0.25* (1-xi (1) *xi (1)) * (1+xi (0)) * (1-xi (2)),
                       0.25* (1-xi (0) *xi (0)) * (1+xi (1)) * (1-xi (2)),
                       0.25* (1-xi (1) *xi (1)) * (1-xi (0)) * (1-xi (2)),
                       0.25* (1-xi (0) *xi (0)) * (1-xi (1)) * (1+xi (2)),
                       0.25* (1-xi (1) *xi (1)) * (1+xi (0)) * (1+xi (2)),
                       0.25* (1-xi (0) *xi (0)) * (1+xi (1)) * (1+xi (2)),
                       0.25* (1-xi (1) *xi (1)) * (1-xi (0)) * (1+xi (2)),
                       0.25* (1-xi (0)) * (1-xi (1)) * (1-xi (2) *xi (2)),
                       0.25* (1+xi (0)) * (1-xi (1)) * (1-xi (2) *xi (2)),
                       0.25* (1+xi (0)) * (1+xi (1)) * (1-xi (2) *xi (2)),
                       0.25* (1-xi (0)) * (1+xi (1)) * (1-xi (2) *xi (2));
                    // clang-format on
                    return N_;
                }
                dNdXiSized dNdXi( const Vector3d& xi ) {
                    dNdXiSized dNdXi_;

                    // clang-format off
                    dNdXi_ <<
                                0.125	*	 (1-xi[1]) *	 (1-xi[2]) *	 (1+2*	xi[0]+xi[1]+xi[2]),
                                -0.125	*	 (1-xi[1]) *	 (1-xi[2]) *	 (1-2*	xi[0]+xi[1]+xi[2]),
                                -0.125	*	 (1+xi[1]) *	 (1-xi[2]) *	 (1-2*	xi[0]-xi[1]+xi[2]),
                                +0.125	*	 (1+xi[1]) *	 (1-xi[2]) *	 (1+2*	xi[0]-xi[1]+xi[2]),
                                +0.125	*	 (1-xi[1]) *	 (1+xi[2]) *	 (1+2*	xi[0]+xi[1]-xi[2]),
                                -0.125	*	 (1-xi[1]) *	 (1+xi[2]) *	 (1-2*	xi[0]+xi[1]-xi[2]),
                                -0.125	*	 (1+xi[1]) *	 (1+xi[2]) *	 (1-2*	xi[0]-xi[1]-xi[2]),
                                +0.125	*	 (1+xi[1]) *	 (1+xi[2]) *	 (1+2*	xi[0]-xi[1]-xi[2]),
                                -0.5	*	xi[0]*	 (1-xi[1]) *	 (1-xi[2]),
                                +0.25	*	 (1-xi[1]*	xi[1]) *	 (1-xi[2]),
                                -0.5	*	xi[0]*	 (1+xi[1]) *	 (1-xi[2]),
                                -0.25	*	 (1-xi[1]*	xi[1]) *	 (1-xi[2]),
                                -0.5	*	xi[0]*	 (1-xi[1]) *	 (1+xi[2]),
                                +0.25	*	 (1-xi[1]*	xi[1]) *	 (1+xi[2]),
                                -0.5	*	xi[0]*	 (1+xi[1]) *	 (1+xi[2]),
                                -0.25	*	 (1-xi[1]*	xi[1]) *	 (1+xi[2]),
                                -0.25	*	 (1-xi[1]) *	 (1-xi[2]*	xi[2]),
                                +0.25	*	 (1-xi[1]) *	 (1-xi[2]*	xi[2]),
                                +0.25	*	 (1+xi[1]) *	 (1-xi[2]*	xi[2]),
                                -0.25	*	 (1+xi[1]) *	 (1-xi[2]*	xi[2]),
                                +0.125	*	 (1-xi[0]) *	 (1-xi[2]) *	 (1+xi[0]+2*	xi[1]+xi[2]),
                                +0.125	*	 (1+xi[0]) *	 (1-xi[2]) *	 (1-xi[0]+2*	xi[1]+xi[2]),
                                -0.125	*	 (1+xi[0]) *	 (1-xi[2]) *	 (1-xi[0]-2*	xi[1]+xi[2]),
                                -0.125	*	 (1-xi[0]) *	 (1-xi[2]) *	 (1+xi[0]-2*	xi[1]+xi[2]),
                                +0.125	*	 (1-xi[0]) *	 (1+xi[2]) *	 (1+xi[0]+2*	xi[1]-xi[2]),
                                +0.125	*	 (1+xi[0]) *	 (1+xi[2]) *	 (1-xi[0]+2*	xi[1]-xi[2]),
                                -0.125	*	 (1+xi[0]) *	 (1+xi[2]) *	 (1-xi[0]-2*	xi[1]-xi[2]),
                                -0.125	*	 (1-xi[0]) *	 (1+xi[2]) *	 (1+xi[0]-2*	xi[1]-xi[2]),
                                -0.25	*	 (1-xi[0]*	xi[0]) *	 (1-xi[2]),
                                -0.5	*	xi[1]*	 (1+xi[0]) *	 (1-xi[2]),
                                +0.25	*	 (1-xi[0]*	xi[0]) *	 (1-xi[2]),
                                -0.5	*	xi[1]*	 (1-xi[0]) *	 (1-xi[2]),
                                -0.25	*	 (1-xi[0]*	xi[0]) *	 (1+xi[2]),
                                -0.5	*	xi[1]*	 (1+xi[0]) *	 (1+xi[2]),
                                +0.25	*	 (1-xi[0]*	xi[0]) *	 (1+xi[2]),
                                -0.5	*	xi[1]*	 (1-xi[0]) *	 (1+xi[2]),
                                -0.25	*	 (1-xi[0]) *	 (1-xi[2]*	xi[2]),
                                -0.25	*	 (1+xi[0]) *	 (1-xi[2]*	xi[2]),
                                +0.25	*	 (1+xi[0]) *	 (1-xi[2]*	xi[2]),
                                +0.25	*	 (1-xi[0]) *	 (1-xi[2]*	xi[2]),
                                +0.125	*	 (1-xi[0]) *	 (1-xi[1]) *	 (1+xi[0]+xi[1]+2*	xi[2]),
                                +0.125	*	 (1+xi[0]) *	 (1-xi[1]) *	 (1-xi[0]+xi[1]+2*	xi[2]),
                                +0.125	*	 (1+xi[0]) *	 (1+xi[1]) *	 (1-xi[0]-xi[1]+2*	xi[2]),
                                +0.125	*	 (1-xi[0]) *	 (1+xi[1]) *	 (1+xi[0]-xi[1]+2*	xi[2]),
                                -0.125	*	 (1-xi[0]) *	 (1-xi[1]) *	 (1+xi[0]+xi[1]-2*	xi[2]),
                                -0.125	*	 (1+xi[0]) *	 (1-xi[1]) *	 (1-xi[0]+xi[1]-2*	xi[2]),
                                -0.125	*	 (1+xi[0]) *	 (1+xi[1]) *	 (1-xi[0]-xi[1]-2*	xi[2]),
                                -0.125	*	 (1-xi[0]) *	 (1+xi[1]) *	 (1+xi[0]-xi[1]-2*	xi[2]),
                                -0.25	*	 (1-xi[0]*	xi[0]) *	 (1-xi[1]),
                                -0.25	*	 (1-xi[1]*	xi[1]) *	 (1+xi[0]),
                                -0.25	*	 (1-xi[0]*	xi[0]) *	 (1+xi[1]),
                                -0.25	*	 (1-xi[1]*	xi[1]) *	 (1-xi[0]),
                                +0.25	*	 (1-xi[0]*	xi[0]) *	 (1-xi[1]),
                                +0.25	*	 (1-xi[1]*	xi[1]) *	 (1+xi[0]),
                                +0.25	*	 (1-xi[0]*	xi[0]) *	 (1+xi[1]),
                                +0.25	*	 (1-xi[1]*	xi[1]) *	 (1-xi[0]),
                                -0.5	*	 (1-xi[0]) *	 (1-xi[1]) *	xi[2],
                                -0.5	*	 (1+xi[0]) *	 (1-xi[1]) *	xi[2],
                                -0.5	*	 (1+xi[0]) *	 (1+xi[1]) *	xi[2],
                                -0.5	*	 (1-xi[0]) *	 (1+xi[1]) *	xi[2];

                    // clang-format on
                    return dNdXi_;
                }

                Matrix3d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates ) {
                    // convenience wrapper to the templated version of the Jacobian
                    return FiniteElement::Jacobian<nDim, nNodes>( dNdXi, coordinates );
                }

                BSized B( const dNdXiSized& dNdX ) {
                    // convenience wrapper to the templated version of the Spatial2D B Operator
                    return FiniteElement::Spatial3D::B<nNodes>( dNdX );
                }

                Vector8i getBoundaryElementIndices( int faceID ) {
                    switch ( faceID ) {
                    case 1: {
                        return ( Vector8i() << 3, 2, 1, 0, 10, 9, 8, 11 ).finished();
                    }
                    case 2: {
                        return ( Vector8i() << 4, 5, 6, 7, 12, 13, 14, 15 ).finished();
                    }
                    case 3: {
                        return ( Vector8i() << 0, 1, 5, 4, 8, 17, 12, 16 ).finished();
                    }
                    case 4: {
                        return ( Vector8i() << 6, 5, 1, 2, 13, 17, 9, 18 ).finished();
                    }
                    case 5: {
                        return ( Vector8i() << 7, 6, 2, 3, 14, 18, 10, 19 ).finished();
                    }
                    case 6: {
                        return ( Vector8i() << 4, 7, 3, 0, 15, 19, 11, 16 ).finished();
                    }
                    default: { throw std::invalid_argument( "Hexa8: invalid face ID specifed" ); }
                    }
                }
            } // end of namespace Hexa20
        }     // end of namespace Spatial3D
    }         // namespace FiniteElement
} // end of namespace bft
