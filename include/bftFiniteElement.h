#pragma once
#include "bftTypedefs.h"
#include <vector>

namespace bft {

    namespace FiniteElement {
        enum ElementShapes {
            Truss2,
            Truss3,
            Quad4,
            Quad8,
            Hexa8,
            Hexa20,
        };

        ElementShapes getElementShapeByMetric( int nDim, int nNodes );

        // 'Expanded' N , aka NBold aka multidimensional Interpolation Operator
        MatrixXd NB( const VectorXd& N,
                     const int       nDoFPerNode ); // Dynamic version

        template <int nDim, int nNodes>
        Matrix<double, nDim, nDim * nNodes> NB( const Matrix<double, 1, nNodes>& N )
        {
            // Alternative Templated version of Interpolation operator NBold;
            Matrix<double, nDim, nDim* nNodes> N_ = Matrix<double, nDim, nDim * nNodes>::Zero();
            for ( int i = 0; i < nNodes; i++ ) {
                for ( int j = 0; j < nDim; j++ ) {
                    N_( j, nDim * i + j ) = N( i );
                }
            }
            return N_;
        }

        MatrixXd Jacobian( const MatrixXd& dN_dXi,
                           const VectorXd& coordinates ); // Dynamic version

        template <int nDim, int nNodes>
        Matrix<double, nDim, nDim> Jacobian( const Matrix<double, nDim, nNodes>&     dNdXi,
                                             const Matrix<double, nDim * nNodes, 1>& coordinates )
        {
            // Alternative Templated version of Jacobian for compile time known
            // sizes
            Matrix<double, nDim, nDim> J_ = Matrix<double, nDim, nDim>::Zero();
            for ( int i = 0; i < nDim; i++ )           // loop over global dimensions
                for ( int j = 0; j < nDim; j++ )       // loop over local dimensions
                    for ( int k = 0; k < nNodes; k++ ) // Loop over nodes
                        J_( i, j ) += dNdXi( j, k ) * coordinates( i + k * nDim );
            return J_;
        }

        VectorXi expandNodeIndicesToCoordinateIndices( const VectorXi& nodeIndices, int nDim );

        namespace Spatial1D {
            namespace Truss2 {

                constexpr int nNodes = 2;
                using NSized         = Matrix<double, 1, nNodes>;
                using dNdXiSized     = Matrix<double, 1, nNodes>;

                NSized     N( double xi );
                dNdXiSized dNdXi( double xi );
            } // namespace Truss2

            namespace Truss3 {

                constexpr int nNodes = 3;
                using NSized         = Matrix<double, 1, nNodes>;
                using dNdXiSized     = Matrix<double, 1, nNodes>;

                NSized     N( double xi );
                dNdXiSized dNdXi( double xi );
            } // namespace Truss3
        }     // namespace Spatial1D

        namespace Spatial2D {
            constexpr int nDim      = 2;
            constexpr int voigtSize = 3;

            template <int nNodes>
            Matrix<double, voigtSize, nNodes * nDim> B( const Matrix<double, nDim, nNodes>& dNdX )
            {

                Matrix<double, voigtSize, nNodes* nDim> B_ = Matrix<double, voigtSize, nNodes * nDim>::Zero();
                for ( int i = 0; i < nNodes; i++ ) {
                    B_( 0, 2 * i )     = dNdX( 0, i );
                    B_( 1, 2 * i + 1 ) = dNdX( 1, i );
                    B_( 2, 2 * i )     = dNdX( 1, i );
                    B_( 2, 2 * i + 1 ) = dNdX( 0, i );
                }
                return B_;
            }

            template <int nNodes>
            Matrix<double, voigtSize, nNodes * nDim> BGreen( const Matrix<double, nDim, nNodes>& dNdX,
                                                             const Matrix2d&                     F )
            {
                // Green-Lagrange Strain Operator for given dNdX and Deformationgradient F
                // Belytschko et. al pp. 213
                Matrix<double, voigtSize, nNodes* nDim> B_ = Matrix<double, voigtSize, nNodes * nDim>::Zero();
                for ( int i = 0; i < nNodes; i++ ) {
                    B_( 0, 2 * i )     = dNdX( 0, i ) * F( 0, 0 );
                    B_( 0, 2 * i + 1 ) = dNdX( 0, i ) * F( 1, 0 );
                    B_( 1, 2 * i )     = dNdX( 1, i ) * F( 0, 1 );
                    B_( 1, 2 * i + 1 ) = dNdX( 1, i ) * F( 1, 1 );
                    B_( 2, 2 * i )     = dNdX( 0, i ) * F( 0, 1 ) + dNdX( 1, i ) * F( 0, 0 );
                    B_( 2, 2 * i + 1 ) = dNdX( 0, i ) * F( 1, 1 ) + dNdX( 1, i ) * F( 1, 0 );
                }
                return B_;
            }

            namespace Quad4 {
                constexpr int nNodes = 4;

                using CoordinateSized = Matrix<double, nNodes * nDim, 1>;
                using NSized          = Matrix<double, 1, nNodes>;
                using dNdXiSized      = Matrix<double, nDim, nNodes>;
                using BSized          = Matrix<double, voigtSize, nNodes * nDim>;

                NSized     N( const Vector2d& xi );
                dNdXiSized dNdXi( const Vector2d& xi );

                // convenience functions; they are wrappers to the corresponding
                // template functions
                Matrix2d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates );
                BSized   B( const dNdXiSized& dNdXi );

                Vector2i getBoundaryElementIndices( int faceID );
            } // namespace Quad4

            namespace Quad8 {
                constexpr int nNodes = 8;

                using CoordinateSized = Matrix<double, nNodes * nDim, 1>;
                using NSized          = Matrix<double, 1, nNodes>;
                using dNdXiSized      = Matrix<double, nDim, nNodes>;
                using BSized          = Matrix<double, voigtSize, nNodes * nDim>;

                NSized     N( const Vector2d& xi );
                dNdXiSized dNdXi( const Vector2d& xi );

                // convenience functions; they are wrappers to the corresponding
                // template functions
                Matrix2d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates );
                BSized   B( const dNdXiSized& dNdXi );

                Vector3i getBoundaryElementIndices( int faceID );
            } // namespace Quad8

        } // end of namespace Spatial2D

        namespace Spatial3D {
            constexpr int nDim      = 3;
            constexpr int voigtSize = 6;

            template <int nNodes>
            Matrix<double, voigtSize, nNodes * nDim> B( const Matrix<double, nDim, nNodes>& dNdX )
            {
                // ABAQUS like notation of strain: ( ε11, ε22, ε33, γ12, γ13, γ23 )
                //   _                                 _
                //  |   dN/dx1    0           0         |
                //  |   0         dN/dx2      0         |
                //  |   0         0           dN/dx3    |
                //  |   dN/dx2    dN/dx1      0         |
                //  |   dN/dx3    0           dN/dx1    |
                //  |_  0         dN/dx3      dN/dx2   _|

                Matrix<double, voigtSize, nNodes* nDim> B_ = Matrix<double, voigtSize, nNodes * nDim>::Zero();

                for ( int i = 0; i < nNodes; i++ ) {
                    B_( 0, nDim * i )     = dNdX( 0, i );
                    B_( 1, nDim * i + 1 ) = dNdX( 1, i );
                    B_( 2, nDim * i + 2 ) = dNdX( 2, i );
                    B_( 3, nDim * i + 0 ) = dNdX( 1, i );
                    B_( 3, nDim * i + 1 ) = dNdX( 0, i );
                    B_( 4, nDim * i + 0 ) = dNdX( 2, i );
                    B_( 4, nDim * i + 2 ) = dNdX( 0, i );
                    B_( 5, nDim * i + 1 ) = dNdX( 2, i );
                    B_( 5, nDim * i + 2 ) = dNdX( 1, i );
                }

                return B_;
            }

            template <int nNodes>
            Matrix<double, voigtSize, nNodes * nDim> BGreen( const Matrix<double, nDim, nNodes>& dNdX,
                                                             const Matrix3d&                     F )
            {
                // Green-Lagrange Strain Operator for given dNdX and Deformationgradient F
                // Belytschko et. al pp. 213

                Matrix<double, voigtSize, nNodes* nDim> B_ = Matrix<double, voigtSize, nNodes * nDim>::Zero();
                for ( int i = 0; i < nNodes; i++ ) {
                    B_( 0, 2 * i )     = dNdX( 0, i ) * F( 0, 0 );
                    B_( 0, 2 * i + 1 ) = dNdX( 0, i ) * F( 1, 0 );
                    B_( 0, 2 * i + 2 ) = dNdX( 0, i ) * F( 2, 0 );

                    B_( 1, 2 * i )     = dNdX( 1, i ) * F( 0, 1 );
                    B_( 1, 2 * i + 1 ) = dNdX( 1, i ) * F( 1, 1 );
                    B_( 1, 2 * i + 2 ) = dNdX( 1, i ) * F( 2, 1 );

                    B_( 2, 2 * i )     = dNdX( 2, i ) * F( 0, 2 );
                    B_( 2, 2 * i + 1 ) = dNdX( 2, i ) * F( 1, 2 );
                    B_( 2, 2 * i + 2 ) = dNdX( 2, i ) * F( 2, 2 );

                    B_( 3, 2 * i )     = dNdX( 0, i ) * F( 0, 1 ) + dNdX( 1, i ) * F( 0, 0 );
                    B_( 3, 2 * i + 1 ) = dNdX( 0, i ) * F( 1, 1 ) + dNdX( 1, i ) * F( 1, 0 );
                    B_( 3, 2 * i + 2 ) = dNdX( 0, i ) * F( 2, 1 ) + dNdX( 1, i ) * F( 2, 0 );

                    B_( 4, 2 * i )     = dNdX( 0, i ) * F( 0, 2 ) + dNdX( 2, i ) * F( 0, 0 );
                    B_( 4, 2 * i + 1 ) = dNdX( 0, i ) * F( 1, 2 ) + dNdX( 2, i ) * F( 1, 0 );
                    B_( 4, 2 * i + 2 ) = dNdX( 0, i ) * F( 2, 2 ) + dNdX( 2, i ) * F( 2, 0 );

                    B_( 5, 2 * i )     = dNdX( 2, i ) * F( 0, 1 ) + dNdX( 1, i ) * F( 0, 2 );
                    B_( 5, 2 * i + 1 ) = dNdX( 2, i ) * F( 1, 1 ) + dNdX( 1, i ) * F( 1, 2 );
                    B_( 5, 2 * i + 2 ) = dNdX( 2, i ) * F( 2, 1 ) + dNdX( 1, i ) * F( 2, 2 );
                }
                return B_;
            }

            namespace Hexa8 {
                constexpr int nNodes  = 8;
                using CoordinateSized = Matrix<double, nNodes * nDim, 1>;
                using NSized          = Matrix<double, 1, nNodes>;
                using dNdXiSized      = Matrix<double, nDim, nNodes>;
                using BSized          = Matrix<double, 6, nNodes * nDim>;

                NSized     N( const Vector3d& xi );
                dNdXiSized dNdXi( const Vector3d& xi );

                // convenience functions; they are wrappers to the corresponding
                // template functions
                Matrix3d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates );
                BSized   B( const dNdXiSized& dNdXi );

                Vector4i getBoundaryElementIndices( int faceID );
            } // namespace Hexa8

            namespace Hexa20 {
                constexpr int nNodes  = 20;
                using CoordinateSized = Matrix<double, nNodes * nDim, 1>;
                using NSized          = Matrix<double, 1, nNodes>;
                using dNdXiSized      = Matrix<double, nDim, nNodes>;
                using BSized          = Matrix<double, 6, nNodes * nDim>;

                NSized     N( const Vector3d& xi );
                dNdXiSized dNdXi( const Vector3d& xi );

                // convenience functions; they are wrappers to the corresponding
                // template functions
                Matrix3d Jacobian( const dNdXiSized& dNdXi, const CoordinateSized& coordinates );
                BSized   B( const dNdXiSized& dNdXi );

                Vector8i getBoundaryElementIndices( int faceID );
            } // namespace Hexa20
        }     // namespace Spatial3D

        class BoundaryElement {

            /* Boundary element, for instance for distributed surface loads
             * */

            struct BoundaryElementGaussPt {
                double   weight;
                VectorXd xi;
                VectorXd N;
                MatrixXd dNdXi;
                MatrixXd J;
                VectorXd normalVector;
                double   integrationArea;
            };

            const int nDim;

            ElementShapes boundaryShape;
            int           nNodes;
            int           nParentCoordinates;

            std::vector<BoundaryElementGaussPt> gaussPts;

            VectorXi boundaryIndicesInParentNodes;
            VectorXi boundaryIndicesInParentCoordinates;
            VectorXd coordinates;

          public:
            BoundaryElement( ElementShapes   parentShape,
                             int             nDim,
                             int             parentFaceNumber,
                             const VectorXd& parentCoordinates );

            VectorXd computeNormalLoadVector();
            VectorXd condenseParentToBoundaryVector( const VectorXd& parentVector );
            VectorXd expandBoundaryToParentVector( const VectorXd& boundaryVector );
        };
    } // namespace FiniteElement

    namespace NumIntegration {
        constexpr double gp2 = 0.577350269189625764509;
        constexpr double gp3 = 0.774596669241483;

        enum IntegrationTypes { FullIntegration, ReducedIntegration };

        struct GaussPtInfo {
            VectorXd xi;
            double   weight;
        };

        const std::vector<GaussPtInfo>& getGaussPointInfo( bft::FiniteElement::ElementShapes shape,
                                                           IntegrationTypes                  integrationType );
        int getNumGaussPoints( bft::FiniteElement::ElementShapes shape, IntegrationTypes integrationType );

        namespace Spatial1D {
            constexpr int nDim = 1;

            // clang-format off
            const std::vector< GaussPtInfo >  gaussPtList1 = {
                { ( VectorXd ( 1 ) << 0 ).finished(),               2.0 }
            };

            const std::vector< GaussPtInfo >  gaussPtList2 = {
                { ( VectorXd ( 1 ) << -gp2 ).finished(),           1.0 },
                { ( VectorXd ( 1 ) << +gp2 ).finished(),           1.0 }
            };

            const std::vector< GaussPtInfo >  gaussPtList3 = {
                { ( VectorXd ( 1 ) << -gp3 ).finished(),            5./9 },
                { ( VectorXd ( 1 ) << 0.   ).finished(),            8./9 },
                { ( VectorXd ( 1 ) << +gp3 ).finished(),            5./9 }
            };
            // clang-format on

        } // namespace Spatial1D

        namespace Spatial2D {
            constexpr int nDim = 2;

            // clang-format off
            const std::vector< GaussPtInfo > gaussPtList1x1 = {
                { Vector2d::Zero(),                             4. }
            };

            const std::vector< GaussPtInfo > gaussPtList2x2 = {
                { ( Vector2d () << +gp2,     +gp2 ).finished(),   1.0 },
                { ( Vector2d () << -gp2,     +gp2 ).finished(),   1.0 },
                { ( Vector2d () << -gp2,     -gp2 ).finished(),   1.0 },
                { ( Vector2d () << +gp2,     -gp2 ).finished(),   1.0 }
            };

            const std::vector< GaussPtInfo > gaussPtList3x3 = {
                { ( Vector2d () << 0,        0.   ).finished(),   64./81},
                { ( Vector2d () << -gp3,     -gp3 ).finished(),   25./81},
                { ( Vector2d () << +gp3,     -gp3 ).finished(),   25./81},
                { ( Vector2d () << +gp3,     +gp3 ).finished(),   25./81},
                { ( Vector2d () << -gp3,     +gp3 ).finished(),   25./81},
                { ( Vector2d () << 0,        -gp3 ).finished(),   40./81},
                { ( Vector2d () << gp3,      0.   ).finished(),   40./81},
                { ( Vector2d () << 0,        +gp3 ).finished(),   40./81},
                { ( Vector2d () << -gp3,     0.   ).finished(),   40./81},
            };
            // clang-format on

            void modifyCharElemLengthAbaqusLike( double& charElemLength, int intPoint );
        } // namespace Spatial2D

        namespace Spatial3D {
            constexpr int nDim = 3;

            // clang-format off
            const std::vector< GaussPtInfo > gaussPtList1x1x1 = {
                { Vector3d::Zero(),                                         8.0 }
            };

            const std::vector< GaussPtInfo > gaussPtList2x2x2 = {
                { ( Vector3d () << -gp2,    -gp2,   -gp2 ).finished(),       1.0},
                { ( Vector3d () << +gp2,    -gp2,   -gp2 ).finished(),       1.0},
                { ( Vector3d () << +gp2,    +gp2,   -gp2 ).finished(),       1.0},
                { ( Vector3d () << -gp2,    +gp2,   -gp2 ).finished(),       1.0},
                { ( Vector3d () << -gp2,    -gp2,   +gp2 ).finished(),       1.0},
                { ( Vector3d () << +gp2,    -gp2,   +gp2 ).finished(),       1.0},
                { ( Vector3d () << +gp2,    +gp2,   +gp2 ).finished(),       1.0},
                { ( Vector3d () << -gp2,    +gp2,   +gp2 ).finished(),       1.0},
            };

            const std::vector< GaussPtInfo > gaussPtList3x3x3 = {
                { ( Vector3d () << -gp3,     -gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << 0,        -gp3,   -gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << +gp3,     -gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << -gp3,     0,      -gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << 0,        0,      -gp3 ).finished(),       0.438957475994513},
                { ( Vector3d () << gp3,      0,      -gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << -gp3,     +gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << 0,        +gp3,   -gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << +gp3,     +gp3,   -gp3 ).finished(),       0.171467764060357},

                { ( Vector3d () << -gp3,     -gp3,   0 ).finished(),          0.274348422496571},
                { ( Vector3d () << 0,        -gp3,   0 ).finished(),          0.438957475994513},
                { ( Vector3d () << +gp3,     -gp3,   0 ).finished(),          0.274348422496571},
                { ( Vector3d () << -gp3,     0,      0 ).finished(),          0.438957475994513},
                { ( Vector3d () << 0,        0,      0 ).finished(),          0.702331961591221},
                { ( Vector3d () << gp3,      0,      0 ).finished(),          0.438957475994513},
                { ( Vector3d () << -gp3,     +gp3,   0 ).finished(),          0.274348422496571},
                { ( Vector3d () << 0,        +gp3,   0 ).finished(),          0.438957475994513},
                { ( Vector3d () << +gp3,     +gp3,   0 ).finished(),          0.274348422496571},

                { ( Vector3d () << -gp3,     -gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << 0,        -gp3,   +gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << +gp3,     -gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << -gp3,     0,      +gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << 0,        0,      +gp3 ).finished(),       0.438957475994513},
                { ( Vector3d () << gp3,      0,      +gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << -gp3,     +gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Vector3d () << 0,        +gp3,   +gp3 ).finished(),       0.274348422496571},
                { ( Vector3d () << +gp3,     +gp3,   +gp3 ).finished(),       0.171467764060357}
            };
            // clang-format on

        } // namespace Spatial3D
    }     // namespace NumIntegration
} // namespace bft
