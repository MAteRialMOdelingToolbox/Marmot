/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Peter Gamnitzer peter.gamnitzer@uibk.ac.at
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotTypedefs.h"
#include <vector>

namespace Marmot {

  namespace FiniteElement {
    enum ElementShapes {
      Bar2,
      Bar3,
      Quad4,
      Quad8,
      Quad9,
      Quad16,
      Tetra4,
      Tetra10,
      Hexa8,
      Hexa20,
      Hexa27,
      Hexa64,
    };

    ElementShapes getElementShapeByMetric( int nDim, int nNodes );

    // 'Expanded' N , aka NBold aka multidimensional Interpolation Operator
    Eigen::MatrixXd NB( const Eigen::VectorXd& N,
                        const int              nDoFPerNode ); // Dynamic version

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nDim, nDim * nNodes > NB( const Eigen::Matrix< double, 1, nNodes >& N )
    {
      // Alternative Templated version of Interpolation operator NBold;
      Eigen::Matrix< double, nDim, nDim* nNodes > N_ = Eigen::Matrix< double, nDim, nDim * nNodes >::Zero();
      for ( int i = 0; i < nNodes; i++ ) {
        for ( int j = 0; j < nDim; j++ ) {
          N_( j, nDim * i + j ) = N( i );
        }
      }
      return N_;
    }

    Eigen::MatrixXd Jacobian( const Eigen::MatrixXd& dN_dXi,
                              const Eigen::VectorXd& coordinates ); // Dynamic version

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nDim, nDim > Jacobian( const Eigen::Matrix< double, nDim, nNodes >&     dNdXi,
                                                  const Eigen::Matrix< double, nDim * nNodes, 1 >& coordinates )
    {
      // Alternative Templated version of Jacobian for compile time known
      // sizes
      Eigen::Matrix< double, nDim, nDim > J_ = Eigen::Matrix< double, nDim, nDim >::Zero();
      for ( int i = 0; i < nDim; i++ )       // loop over global dimensions
        for ( int j = 0; j < nDim; j++ )     // loop over local dimensions
          for ( int k = 0; k < nNodes; k++ ) // Loop over nodes
            J_( i, j ) += dNdXi( j, k ) * coordinates( i + k * nDim );
      return J_;
    }

    Eigen::VectorXi expandNodeIndicesToCoordinateIndices( const Eigen::VectorXi& nodeIndices, int nDim );

    namespace Spatial1D {
      namespace Bar2 {

        constexpr int nNodes = 2;
        using NSized         = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized     = Eigen::Matrix< double, 1, nNodes >;

        NSized     N( double xi );
        dNdXiSized dNdXi( double xi );
      } // namespace Bar2

      namespace Bar3 {

        constexpr int nNodes = 3;
        using NSized         = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized     = Eigen::Matrix< double, 1, nNodes >;

        NSized     N( double xi );
        dNdXiSized dNdXi( double xi );
      } // namespace Bar3
    }   // namespace Spatial1D

    namespace Spatial2D {
      constexpr int nDim      = 2;
      constexpr int voigtSize = 3;

      template < int nNodes >
      Eigen::Matrix< double, voigtSize, nNodes * nDim > B( const Eigen::Matrix< double, nDim, nNodes >& dNdX )
      {

        Eigen::Matrix< double, voigtSize, nNodes* nDim > B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();
        for ( int i = 0; i < nNodes; i++ ) {
          B_( 0, 2 * i )     = dNdX( 0, i );
          B_( 1, 2 * i + 1 ) = dNdX( 1, i );
          B_( 2, 2 * i )     = dNdX( 1, i );
          B_( 2, 2 * i + 1 ) = dNdX( 0, i );
        }
        return B_;
      }
      namespace axisymmetric {
        constexpr int voigtSize = 4;

        template < int nNodes >
        Eigen::Matrix< double, voigtSize, nNodes * nDim > B( const Eigen::Matrix< double, nDim, nNodes >& dNdX,
                                                             const Eigen::Matrix< double, 1, nNodes >&    N,
                                                             const Eigen::Matrix< double, nDim, 1 >&      x_gauss )
        {

          Eigen::Matrix< double, voigtSize, nNodes* nDim >
            B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();
          for ( int i = 0; i < nNodes; i++ ) {
            B_( 0, 2 * i )     = dNdX( 0, i );
            B_( 1, 2 * i + 1 ) = dNdX( 1, i );
            B_( 2, 2 * i )     = N( i ) / x_gauss( 0 ); // ( N_I / r )
            B_( 3, 2 * i )     = dNdX( 1, i );
            B_( 3, 2 * i + 1 ) = dNdX( 0, i );
          }
          return B_;
        }
      } // namespace axisymmetric

      template < int nNodes >
      Eigen::Matrix< double, voigtSize, nNodes * nDim > BGreen( const Eigen::Matrix< double, nDim, nNodes >& dNdX,
                                                                const Eigen::Matrix2d&                       F )
      {
        // Green-Lagrange Strain Operator for given dNdX and Deformationgradient F
        // Belytschko et. al pp. 213
        Eigen::Matrix< double, voigtSize, nNodes* nDim > B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();
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

        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector2d& xi );
        dNdXiSized dNdXi( const Eigen::Vector2d& xi );

        Eigen::Vector2i getBoundaryElementIndices( int faceID );
      } // namespace Quad4

      namespace Quad8 {
        constexpr int nNodes = 8;

        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector2d& xi );
        dNdXiSized dNdXi( const Eigen::Vector2d& xi );

        Eigen::Vector3i getBoundaryElementIndices( int faceID );
      } // namespace Quad8

    }   // end of namespace Spatial2D

    namespace Spatial3D {
      constexpr int nDim      = 3;
      constexpr int voigtSize = 6;

      template < int nNodes >
      Eigen::Matrix< double, voigtSize, nNodes * nDim > B( const Eigen::Matrix< double, nDim, nNodes >& dNdX )
      {
        // ABAQUS like notation of strain: ( ε11, ε22, ε33, γ12, γ13, γ23 )
        //   _                                 _
        //  |   dN/dx1    0           0         |
        //  |   0         dN/dx2      0         |
        //  |   0         0           dN/dx3    |
        //  |   dN/dx2    dN/dx1      0         |
        //  |   dN/dx3    0           dN/dx1    |
        //  |_  0         dN/dx3      dN/dx2   _|

        Eigen::Matrix< double, voigtSize, nNodes* nDim > B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();

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

      template < int nNodes >
      Eigen::Matrix< double, voigtSize, nNodes * nDim > B_bar( const Eigen::Matrix< double, nDim, nNodes >& dNdX,
                                                               const Eigen::Matrix< double, nDim, nNodes >& dNdX0 )
      {
        // ABAQUS like notation of strain: ( ε11, ε22, ε33, γ12, γ13, γ23 )

        // Implementation of selective reduced integration using the B-bar method (T.J.R. Hughes, 1980).
        // The matrix dNdX0 is evaluated at the center of the element (Xi_i = 0).
        // The B-bar method modifies the strain-displacement matrix (B matrix) such that the volumetric part
        // of the strain is replaced by its average over all quadrature points. This technique helps to
        // alleviate volumetric locking in nearly incompressible materials.
        // The modified B matrix for the node a is defined as:
        /*
        B̄_a =
            [ B5   B6   B8 ]
            [ B4   B7   B8 ]
            [ B4   B6   B9 ]
            [ B2   B1    0 ]
            [  0   B3   B2 ]
            [ B3    0   B1 ]
        */

        Eigen::Matrix< double, voigtSize, nNodes* nDim > B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();

        double B1, B2, B3, B4, B5, B6, B7, B8, B9;
        double B1_bar, B2_bar, B3_bar;

        for ( int i = 0; i < nNodes; i++ ) {
          B1                    = dNdX( 0, i );
          B2                    = dNdX( 1, i );
          B3                    = dNdX( 2, i );
          B1_bar                = dNdX0( 0, i );
          B2_bar                = dNdX0( 1, i );
          B3_bar                = dNdX0( 2, i );
          B4                    = ( B1_bar - B1 ) / 3.;
          B5                    = B1 + B4;
          B6                    = ( B2_bar - B2 ) / 3.;
          B7                    = B2 + B6;
          B8                    = ( B3_bar - B3 ) / 3.;
          B9                    = B3 + B8;
          B_( 0, nDim * i )     = B5;
          B_( 0, nDim * i + 1 ) = B6;
          B_( 0, nDim * i + 2 ) = B8;
          B_( 1, nDim * i )     = B4;
          B_( 1, nDim * i + 1 ) = B7;
          B_( 1, nDim * i + 2 ) = B8;
          B_( 2, nDim * i )     = B4;
          B_( 2, nDim * i + 1 ) = B6;
          B_( 2, nDim * i + 2 ) = B9;

          // shear part is the same as in the normal B matrix
          B_( 3, nDim * i + 0 ) = dNdX( 1, i );
          B_( 3, nDim * i + 1 ) = dNdX( 0, i );
          B_( 4, nDim * i + 0 ) = dNdX( 2, i );
          B_( 4, nDim * i + 2 ) = dNdX( 0, i );
          B_( 5, nDim * i + 1 ) = dNdX( 2, i );
          B_( 5, nDim * i + 2 ) = dNdX( 1, i );
        }

        return B_;
      }

      template < int nNodes >
      Eigen::Matrix< double, voigtSize, nNodes * nDim > BGreen( const Eigen::Matrix< double, nDim, nNodes >& dNdX,
                                                                const Eigen::Matrix3d&                       F )
      {
        // Green-Lagrange Strain Operator for given dNdX and Deformationgradient F
        // Belytschko et. al pp. 213

        Eigen::Matrix< double, voigtSize, nNodes* nDim > B_ = Eigen::Matrix< double, voigtSize, nNodes * nDim >::Zero();
        for ( int i = 0; i < nNodes; i++ ) {
          B_( 0, nDim * i )     = dNdX( 0, i ) * F( 0, 0 );
          B_( 0, nDim * i + 1 ) = dNdX( 0, i ) * F( 1, 0 );
          B_( 0, nDim * i + 2 ) = dNdX( 0, i ) * F( 2, 0 );

          B_( 1, nDim * i )     = dNdX( 1, i ) * F( 0, 1 );
          B_( 1, nDim * i + 1 ) = dNdX( 1, i ) * F( 1, 1 );
          B_( 1, nDim * i + 2 ) = dNdX( 1, i ) * F( 2, 1 );

          B_( 2, nDim * i )     = dNdX( 2, i ) * F( 0, 2 );
          B_( 2, nDim * i + 1 ) = dNdX( 2, i ) * F( 1, 2 );
          B_( 2, nDim * i + 2 ) = dNdX( 2, i ) * F( 2, 2 );

          B_( 3, nDim * i )     = dNdX( 0, i ) * F( 0, 1 ) + dNdX( 1, i ) * F( 0, 0 );
          B_( 3, nDim * i + 1 ) = dNdX( 0, i ) * F( 1, 1 ) + dNdX( 1, i ) * F( 1, 0 );
          B_( 3, nDim * i + 2 ) = dNdX( 0, i ) * F( 2, 1 ) + dNdX( 1, i ) * F( 2, 0 );

          B_( 4, nDim * i )     = dNdX( 0, i ) * F( 0, 2 ) + dNdX( 2, i ) * F( 0, 0 );
          B_( 4, nDim * i + 1 ) = dNdX( 0, i ) * F( 1, 2 ) + dNdX( 2, i ) * F( 1, 0 );
          B_( 4, nDim * i + 2 ) = dNdX( 0, i ) * F( 2, 2 ) + dNdX( 2, i ) * F( 2, 0 );

          B_( 5, nDim * i )     = dNdX( 2, i ) * F( 0, 1 ) + dNdX( 1, i ) * F( 0, 2 );
          B_( 5, nDim * i + 1 ) = dNdX( 2, i ) * F( 1, 1 ) + dNdX( 1, i ) * F( 1, 2 );
          B_( 5, nDim * i + 2 ) = dNdX( 2, i ) * F( 2, 1 ) + dNdX( 1, i ) * F( 2, 2 );
        }
        return B_;
      }

      namespace Tetra4 {

        constexpr int nNodes  = 4;
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector3d& xi );
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        Eigen::Vector3i getBoundaryElementIndices( int faceID );

      } // namespace Tetra4

      namespace Tetra10 {

        constexpr int nNodes  = 10;
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector3d& xi );
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        Eigen::Vector3i getBoundaryElementIndices( int faceID );

      } // namespace Tetra10

      namespace Hexa8 {
        constexpr int nNodes  = 8;
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector3d& xi );
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        Eigen::Vector4i getBoundaryElementIndices( int faceID );
      } // namespace Hexa8

      namespace Hexa20 {
        constexpr int nNodes  = 20;
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >;
        using NSized          = Eigen::Matrix< double, 1, nNodes >;
        using dNdXiSized      = Eigen::Matrix< double, nDim, nNodes >;

        NSized     N( const Eigen::Vector3d& xi );
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        Marmot::Vector8i getBoundaryElementIndices( int faceID );
      } // namespace Hexa20
    }   // namespace Spatial3D

    class BoundaryElement {

      /* Boundary element, for instance for distributed surface loads
       * */

      struct BoundaryElementQuadraturePoint {
        double          weight;
        double          JxW;
        Eigen::VectorXd xi;
        Eigen::VectorXd N;
        Eigen::MatrixXd dNdXi;
        Eigen::MatrixXd dx_dXi;
        Eigen::VectorXd areaVector;
      };

      const int nDim;

      ElementShapes boundaryShape;
      int           nNodes;
      int           nParentCoordinates;

      std::vector< BoundaryElementQuadraturePoint > quadraturePoints;

      Eigen::VectorXi mapBoundaryToParentScalar;
      Eigen::VectorXi mapBoundaryToParentVectorial;
      Eigen::VectorXd coordinates;

    public:
      BoundaryElement( ElementShapes          parentShape,
                       int                    nDim,
                       int                    parentFaceNumber,
                       const Eigen::VectorXd& parentCoordinates );

      Eigen::VectorXd computeScalarLoadVector();
      Eigen::MatrixXd computeDScalarLoadVector_dCoordinates();

      /// compute the element load vector for a unit vectorial load normal to the surface.
      Eigen::VectorXd computeSurfaceNormalVectorialLoadVector();
      Eigen::VectorXd computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements();
      Eigen::MatrixXd computeDSurfaceNormalVectorialLoadVector_dCoordinates();

      /// compute the element load vector for a unit vectorial load in a given direction.
      Eigen::VectorXd computeVectorialLoadVector( const Eigen::VectorXd& direction );
      Eigen::MatrixXd computeDVectorialLoadVector_dCoordinates( const Eigen::VectorXd& direction );

      Eigen::VectorXd condenseParentToBoundaryScalar( const Eigen::VectorXd& parentVector );
      void            assembleIntoParentScalar( const Eigen::VectorXd&        boundaryVector,
                                                Eigen::Ref< Eigen::VectorXd > ParentVector );
      void assembleIntoParentStiffnessScalar( const Eigen::MatrixXd& KBoundary, Eigen::Ref< Eigen::MatrixXd > KParent );

      Eigen::VectorXd condenseParentToBoundaryVectorial( const Eigen::VectorXd& parentVector );
      void            assembleIntoParentVectorial( const Eigen::VectorXd&        boundaryVector,
                                                   Eigen::Ref< Eigen::VectorXd > ParentVector );
      void            assembleIntoParentStiffnessVectorial( const Eigen::MatrixXd&        KBoundary,
                                                            Eigen::Ref< Eigen::MatrixXd > KParent );
    };
  } // namespace FiniteElement

  namespace FiniteElement::Quadrature {
    constexpr double gp2 = 0.577350269189625764509;
    constexpr double gp3 = 0.774596669241483;

    enum IntegrationTypes { FullIntegration, ReducedIntegration };

    struct QuadraturePointInfo {
      Eigen::VectorXd xi;
      double          weight;
    };

    const std::vector< QuadraturePointInfo >& getGaussPointInfo( Marmot::FiniteElement::ElementShapes shape,
                                                                 IntegrationTypes                     integrationType );

    int getNumGaussPoints( Marmot::FiniteElement::ElementShapes shape, IntegrationTypes integrationType );

    namespace Spatial1D {
      constexpr int nDim = 1;

      // clang-format off
            const std::vector< QuadraturePointInfo >  gaussPointList1 = {
                { ( Eigen::VectorXd ( 1 ) << 0 ).finished(),               2.0 }
            };

            const std::vector< QuadraturePointInfo >  gaussPointList2 = {
                { ( Eigen::VectorXd ( 1 ) << -gp2 ).finished(),           1.0 },
                { ( Eigen::VectorXd ( 1 ) << +gp2 ).finished(),           1.0 }
            };

            const std::vector< QuadraturePointInfo >  gaussPointList3 = {
                { ( Eigen::VectorXd ( 1 ) << -gp3 ).finished(),            5./9 },
                { ( Eigen::VectorXd ( 1 ) << 0.   ).finished(),            8./9 },
                { ( Eigen::VectorXd ( 1 ) << +gp3 ).finished(),            5./9 }
            };
      // clang-format on

    } // namespace Spatial1D

    namespace Spatial2D {
      constexpr int nDim = 2;

      // clang-format off
            const std::vector< QuadraturePointInfo > gaussPointList1x1 = {
                { Eigen::Vector2d::Zero(),                             4. }
            };

            const std::vector< QuadraturePointInfo > gaussPointList2x2 = {
                { ( Eigen::Vector2d () << +gp2,     +gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << -gp2,     +gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << -gp2,     -gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << +gp2,     -gp2 ).finished(),   1.0 }
            };

            const std::vector< QuadraturePointInfo > gaussPointList3x3 = {
                { ( Eigen::Vector2d () << 0,        0.   ).finished(),   64./81},
                { ( Eigen::Vector2d () << -gp3,     -gp3 ).finished(),   25./81},
                { ( Eigen::Vector2d () << +gp3,     -gp3 ).finished(),   25./81},
                { ( Eigen::Vector2d () << +gp3,     +gp3 ).finished(),   25./81},
                { ( Eigen::Vector2d () << -gp3,     +gp3 ).finished(),   25./81},
                { ( Eigen::Vector2d () << 0,        -gp3 ).finished(),   40./81},
                { ( Eigen::Vector2d () << gp3,      0.   ).finished(),   40./81},
                { ( Eigen::Vector2d () << 0,        +gp3 ).finished(),   40./81},
                { ( Eigen::Vector2d () << -gp3,     0.   ).finished(),   40./81},
            };
      // clang-format on

      void modifyCharElemLengthAbaqusLike( double& charElemLength, int intPoint );
    } // namespace Spatial2D

    namespace Spatial3D {
      constexpr int nDim = 3;

      // clang-format off
            const inline std::vector< QuadraturePointInfo > gaussPointList1x1x1 = {
                { Eigen::Vector3d::Zero(),                                         8.0 }
            };

            const inline std::vector< QuadraturePointInfo > gaussPointListTetra4 = {
                { (Eigen::Vector3d() << 1./4, 1./4, 1./4).finished(),  1./6}
            };

            const inline std::vector< QuadraturePointInfo > gaussPointListTetra10 = {
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20     ).finished(),  1./24},
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20,    (5+3*std::sqrt(5))/20   ).finished(),  1./24},
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5+3*std::sqrt(5))/20,  (5-std::sqrt(5))/20     ).finished(),  1./24},
                { (Eigen::Vector3d() << (5+3*std::sqrt(5))/20,  (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20     ).finished(),  1./24},
            };

            const inline std::vector< QuadraturePointInfo > gaussPointList2x2x2 = {
                { ( Eigen::Vector3d () << -gp2,    -gp2,   -gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << +gp2,    -gp2,   -gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << +gp2,    +gp2,   -gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << -gp2,    +gp2,   -gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << -gp2,    -gp2,   +gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << +gp2,    -gp2,   +gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << +gp2,    +gp2,   +gp2 ).finished(),       1.0},
                { ( Eigen::Vector3d () << -gp2,    +gp2,   +gp2 ).finished(),       1.0},
            };

            const inline std::vector< QuadraturePointInfo > gaussPointList3x3x3 = {
                { ( Eigen::Vector3d () << -gp3,     -gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << 0,        -gp3,   -gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << +gp3,     -gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << -gp3,     0,      -gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << 0,        0,      -gp3 ).finished(),       0.438957475994513},
                { ( Eigen::Vector3d () << gp3,      0,      -gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << -gp3,     +gp3,   -gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << 0,        +gp3,   -gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << +gp3,     +gp3,   -gp3 ).finished(),       0.171467764060357},

                { ( Eigen::Vector3d () << -gp3,     -gp3,   0 ).finished(),          0.274348422496571},
                { ( Eigen::Vector3d () << 0,        -gp3,   0 ).finished(),          0.438957475994513},
                { ( Eigen::Vector3d () << +gp3,     -gp3,   0 ).finished(),          0.274348422496571},
                { ( Eigen::Vector3d () << -gp3,     0,      0 ).finished(),          0.438957475994513},
                { ( Eigen::Vector3d () << 0,        0,      0 ).finished(),          0.702331961591221},
                { ( Eigen::Vector3d () << gp3,      0,      0 ).finished(),          0.438957475994513},
                { ( Eigen::Vector3d () << -gp3,     +gp3,   0 ).finished(),          0.274348422496571},
                { ( Eigen::Vector3d () << 0,        +gp3,   0 ).finished(),          0.438957475994513},
                { ( Eigen::Vector3d () << +gp3,     +gp3,   0 ).finished(),          0.274348422496571},

                { ( Eigen::Vector3d () << -gp3,     -gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << 0,        -gp3,   +gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << +gp3,     -gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << -gp3,     0,      +gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << 0,        0,      +gp3 ).finished(),       0.438957475994513},
                { ( Eigen::Vector3d () << gp3,      0,      +gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << -gp3,     +gp3,   +gp3 ).finished(),       0.171467764060357},
                { ( Eigen::Vector3d () << 0,        +gp3,   +gp3 ).finished(),       0.274348422496571},
                { ( Eigen::Vector3d () << +gp3,     +gp3,   +gp3 ).finished(),       0.171467764060357}
            };
      // clang-format on

    } // namespace Spatial3D
  }   // namespace FiniteElement::Quadrature
} // namespace Marmot
