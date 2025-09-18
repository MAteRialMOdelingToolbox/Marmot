#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>

using namespace Eigen;
namespace Marmot {
  namespace FiniteElement {
    namespace Spatial3D {
      namespace Tetra4 {

        NSized N( const Vector3d& xi )
        {
          /*
           *            xi(1)
           *
           *            (3)
           *            /|  \
           *           / |    \
           *          /  |      \
           *         /   |        \
           *        /   (1)_________\(2)   xi(0)
           *       /    ´        _-
           *      /   ´     _-
           *     /  ´   _-
           *    (4)_-
           *
           *  xi(2)
           *  */

          NSized N_;

          // clang-format off
                    N_ << 1 - xi(0) - xi(1) - xi(2),
                       xi(0),
                       xi(1),
                       xi(2);

          // clang-format on

          return N_;
        }
        dNdXiSized dNdXi( const Vector3d& xi )
        {
          dNdXiSized dNdXi_;

          // clang-format off
                    dNdXi_ << -1,
                           1,
                           0,
                           0,

                           -1,
                           0,
                           1,
                           0,

                           -1,
                           0,
                           0,
                           1;

          // clang-format on

          return dNdXi_;
        }

        Vector3i getBoundaryElementIndices( int faceID )
        {
          throw std::invalid_argument( "Tetra4: invalid face ID specifed" );
        }

      } // namespace Tetra4
      namespace Tetra10 {

        NSized N( const Vector3d& xi )
        {
          /*
           *            xi(1)
           *
           *            (2)
           *            /|  \
           *           / |    \
           *          / (6)   (5)
           *       (9)   |        \
           *        /   (3)__(7)____\(1)   xi(0)
           *       /    ´        _-
           *      /  (10)   (8)
           *     /  ´   _-
           *    (4)_-
           *
           *  xi(2)
           *
           *  */

          NSized N_;

          // clang-format off
                    N_ <<   xi( 1 ) * ( 2 * xi( 1 ) - 1 ),
                            xi( 2 ) * ( 2 * xi( 2 ) - 1 ),
                            ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ) * ( 1 - 2 * xi( 0 ) - 2 * xi( 1 ) - 2 * xi( 2 ) ),
                            xi( 0 ) * ( 2 * xi( 0 ) - 1 ),
                            4 * xi( 1 ) * xi( 2 ),
                            4 * xi( 2 ) * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ),
                            4 * xi( 1 ) * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ),
                            4 * xi( 0 ) * xi( 1 ),
                            4 * xi( 0 ) * xi( 2 ),
                            4 * xi( 0 ) * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) );
          // clang-format on

          return N_;
        }
        dNdXiSized dNdXi( const Vector3d& xi )
        {
          dNdXiSized dNdXi_;

          // clang-format off
                    dNdXi_ <<
                            0,
                            0,
                            ( -1 ) * ( 1 - 2 * xi( 0 ) - 2 * xi( 1 ) - 2 * xi( 2 ) )
                                +
                            ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ) * (-2 ),
                            1 * ( 2 * xi( 0 ) - 1 )
                                +
                            xi( 0 ) * ( 2 ),
                            0,
                            4 * xi( 2 ) * ( -1 ),
                            4 * xi( 1 ) * ( -1 ),
                            4 * 1 * xi( 1 ),
                            4 * 1 * xi( 2 ),
                            4 * 1 * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) )
                                +
                                4 * xi( 0 ) * ( -1 ),


                            1 * ( 2 * xi( 1 ) - 1 )
                                +
                            xi( 1 ) * ( 2 ),
                            0,
                            ( -1 ) * ( 1 - 2 * xi( 0 ) - 2 * xi( 1 ) - 2 * xi( 2 ) )
                                +
                            ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ) * ( -2 ),
                            0,
                            4 * 1 * xi( 2 ),
                            4 * xi( 2 ) * ( -1 ),
                            4 * 1 * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) )
                                +
                            4 * xi( 1 ) * ( -1 ),
                            4 * xi( 0 ) * 1,
                            0,
                            4 * xi( 0 ) * ( -1 ),


                            0,
                            1 * ( 2 * xi( 2 ) - 1 )
                                +
                            xi( 2 ) * ( 2 ),
                            ( -1 ) * ( 1 - 2 * xi( 0 ) - 2 * xi( 1 ) - 2 * xi( 2 ) )
                                +
                            ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ) * ( -2 ),
                            0,
                            4 * xi( 1 ) * 1,
                            4 * xi( 2 ) * ( -1 )
                                +
                            4 * 1 * ( 1 - xi( 0 ) - xi( 1 ) - xi( 2 ) ),
                            4 * xi( 1 ) * ( -1 ),
                            0,
                            4 * xi( 0 ) * 1,
                            4 * xi( 0 ) * ( -1 );

          // clang-format on

          return dNdXi_;
        }

        Vector3i getBoundaryElementIndices( int faceID )
        {
          throw std::invalid_argument( "Tetra4: invalid face ID specifed" );
        }

      } // namespace Tetra10

      namespace Hexa8 {

        NSized N( const Vector3d& xi )
        {
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

          N_ *= 1. / 8;
          return N_;
        }
        dNdXiSized dNdXi( const Vector3d& xi )
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

        Vector4i getBoundaryElementIndices( int faceID )
        {
          switch ( faceID ) {
          // bottom
          case 1: {
            return ( Vector4i() << 3, 2, 1, 0 ).finished();
          }
            // top
          case 2: {
            return ( Vector4i() << 4, 5, 6, 7 ).finished();
          }
            // front
          case 3: {
            return ( Vector4i() << 0, 1, 5, 4 ).finished();
          }
            // right
          case 4: {
            return ( Vector4i() << 6, 5, 1, 2 ).finished();
          }
            // back
          case 5: {
            return ( Vector4i() << 7, 6, 2, 3 ).finished();
          }
            // left
          case 6: {
            return ( Vector4i() << 4, 7, 3, 0 ).finished();
          }
          default: {
            throw std::invalid_argument( "Hexa8: invalid face ID specifed" );
          }
          }
        }

      } // end of namespace Hexa8

      namespace Hexa20 {
        NSized N( const Vector3d& xi )
        {
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
        dNdXiSized dNdXi( const Vector3d& xi )
        {
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

        Vector8i getBoundaryElementIndices( int faceID )
        {
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
          default: {
            throw std::invalid_argument( "Hexa8: invalid face ID specifed" );
          }
          }
        }
      } // end of namespace Hexa20
    }   // end of namespace Spatial3D
  }     // namespace FiniteElement
} // end of namespace Marmot
