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

/**
 * @file MarmotFiniteElement.h
 * @brief Core finite element utilities: element shapes, shape functions, Jacobian,
 *        and B-operator for 1-D, 2-D and 3-D isoparametric elements.
 */

namespace Marmot {

  namespace FiniteElement {
    /**
     * @brief Supported isoparametric element shapes.
     */
    enum ElementShapes {
      Bar2,    ///< 1-D 2-node bar
      Bar3,    ///< 1-D 3-node bar
      Quad4,   ///< 2-D 4-node quadrilateral
      Quad8,   ///< 2-D 8-node serendipity quadrilateral
      Quad9,   ///< 2-D 9-node quadrilateral
      Quad16,  ///< 2-D 16-node quadrilateral
      Tetra4,  ///< 3-D 4-node tetrahedron
      Tetra10, ///< 3-D 10-node tetrahedron
      Hexa8,   ///< 3-D 8-node hexahedron
      Hexa20,  ///< 3-D 20-node serendipity hexahedron
      Hexa27,  ///< 3-D 27-node hexahedron
      Hexa64,  ///< 3-D 64-node hexahedron
    };

    /**
     * @brief Determine the element shape from spatial dimension and node count.
     * @param[in] nDim   Number of spatial dimensions.
     * @param[in] nNodes Number of element nodes.
     * @return Matching @ref ElementShapes enumerator.
     */
    ElementShapes getElementShapeByMetric( int nDim, int nNodes );

    /**
     * @brief Compute the expanded interpolation operator N_B (dynamic version).
     * @param[in] N           Shape function row vector.
     * @param[in] nDoFPerNode Number of degrees of freedom per node.
     * @return Expanded matrix of size @c nDoFPerNode × (nNodes * nDoFPerNode).
     */
    Eigen::MatrixXd NB( const Eigen::VectorXd& N,
                        const int              nDoFPerNode ); // Dynamic version

    /**
     * @brief Compute the expanded interpolation operator N_B (compile-time size version).
     * @tparam nDim   Number of spatial dimensions.
     * @tparam nNodes Number of element nodes.
     * @param[in] N Shape function row vector.
     * @return Expanded matrix of size @c nDim × (nDim * nNodes).
     */
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

    /**
     * @brief Compute the element Jacobian matrix (dynamic version).
     * @param[in] dN_dXi    Matrix of shape function natural derivatives.
     * @param[in] coordinates Flat vector of nodal coordinates.
     * @return Jacobian matrix.
     */
    Eigen::MatrixXd Jacobian( const Eigen::MatrixXd& dN_dXi,
                              const Eigen::VectorXd& coordinates ); // Dynamic version

    /**
     * @brief Compute the element Jacobian matrix (compile-time size version).
     * @tparam nDim   Number of spatial dimensions.
     * @tparam nNodes Number of element nodes.
     * @param[in] dNdXi       Matrix of shape function natural derivatives.
     * @param[in] coordinates Flat nodal-coordinate vector.
     * @return Square Jacobian matrix of size @c nDim × nDim.
     */
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

    /**
     * @brief Expand scalar node indices to coordinate (DOF) indices.
     * @param[in] nodeIndices Vector of node indices.
     * @param[in] nDim        Number of spatial dimensions (DOFs per node).
     * @return Vector of coordinate indices.
     */
    Eigen::VectorXi expandNodeIndicesToCoordinateIndices( const Eigen::VectorXi& nodeIndices, int nDim );

    namespace Spatial1D {
      namespace Bar2 {

        constexpr int nNodes = 2;                                  ///< Number of nodes for a Bar2 element.
        using NSized         = Eigen::Matrix< double, 1, nNodes >; ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, 1, nNodes >; ///< Row vector of shape function natural derivatives.

        /// @brief Evaluate Bar2 shape functions at natural coordinate @p xi.
        NSized N( double xi );
        /// @brief Evaluate Bar2 shape function derivatives at natural coordinate @p xi.
        dNdXiSized dNdXi( double xi );
      } // namespace Bar2

      namespace Bar3 {

        constexpr int nNodes = 3;                                  ///< Number of nodes for a Bar3 element.
        using NSized         = Eigen::Matrix< double, 1, nNodes >; ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, 1, nNodes >; ///< Row vector of shape function natural derivatives.

        /// @brief Evaluate Bar3 shape functions at natural coordinate @p xi.
        NSized N( double xi );
        /// @brief Evaluate Bar3 shape function derivatives at natural coordinate @p xi.
        dNdXiSized dNdXi( double xi );
      } // namespace Bar3
    }   // namespace Spatial1D

    namespace Spatial2D {
      constexpr int nDim      = 2; ///< Number of spatial dimensions for 2-D elements.
      constexpr int voigtSize = 3; ///< Voigt vector size for 2-D plane elements (σ₁₁, σ₂₂, σ₁₂).

      /**
       * @brief Compute the 2-D linear B-operator from physical derivatives.
       * @tparam nNodes Number of element nodes.
       * @param[in] dNdX Physical shape function derivatives (2 × nNodes).
       * @return B-operator matrix of size 3 × (2 * nNodes).
       */
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
      namespace Axisymmetric {
        constexpr int voigtSize = 4; ///< Voigt vector size for axisymmetric elements (σ_rr, σ_zz, σ_θθ, σ_rz).

        /**
         * @brief Compute the axisymmetric B-operator from physical derivatives.
         * @tparam nNodes Number of element nodes.
         * @param[in] dNdX    Physical shape function derivatives (2 × nNodes).
         * @param[in] N       Shape function values (1 × nNodes).
         * @param[in] x_gauss Physical coordinates of the Gauss point (r, z).
         * @return B-operator matrix of size 4 × (2 * nNodes).
         */
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
      } // namespace Axisymmetric

      /**
       * @brief Compute the 2-D Green–Lagrange B-operator.
       * @tparam nNodes Number of element nodes.
       * @param[in] dNdX Physical shape function derivatives (2 × nNodes).
       * @param[in] F    2×2 deformation gradient.
       * @return Green–Lagrange B-operator of size 3 × (2 * nNodes).
       */
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
        constexpr int nNodes = 4;                                          ///< Number of nodes for a Quad4 element.

        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Quad4 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector2d& xi );
        /// @brief Evaluate Quad4 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector2d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Eigen::Vector2i getBoundaryElementIndices( int faceID );
      } // namespace Quad4

      namespace Quad8 {
        constexpr int nNodes = 8;                                          ///< Number of nodes for a Quad8 element.

        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Quad8 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector2d& xi );
        /// @brief Evaluate Quad8 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector2d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Eigen::Vector3i getBoundaryElementIndices( int faceID );
      } // namespace Quad8

    }   // namespace Spatial2D

    namespace Spatial3D {
      constexpr int nDim = 3; ///< Number of spatial dimensions for 3-D elements.
      constexpr int voigtSize = 6; ///< Voigt vector size for 3-D elements (σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃).

      /**
       * @brief Compute the 3-D linear B-operator from physical derivatives.
       * @tparam nNodes Number of element nodes.
       * @param[in] dNdX Physical shape function derivatives (3 × nNodes).
       * @return B-operator matrix of size 6 × (3 * nNodes).
       */
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

      /**
       * @brief Compute the 3-D B̄-operator using the B-bar method (Hughes, 1980).
       * @tparam nNodes Number of element nodes.
       * @param[in] dNdX  Physical shape function derivatives at the integration point.
       * @param[in] dNdX0 Physical shape function derivatives at the element centre.
       * @return Modified B-operator matrix of size 6 × (3 * nNodes).
       */
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

      /**
       * @brief Compute the 3-D Green–Lagrange B-operator.
       * @tparam nNodes Number of element nodes.
       * @param[in] dNdX Physical shape function derivatives (3 × nNodes).
       * @param[in] F    3×3 deformation gradient.
       * @return Green–Lagrange B-operator of size 6 × (3 * nNodes).
       */
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

        constexpr int nNodes  = 4;                                         ///< Number of nodes for a Tetra4 element.
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Tetra4 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector3d& xi );
        /// @brief Evaluate Tetra4 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Eigen::Vector3i getBoundaryElementIndices( int faceID );

      } // namespace Tetra4

      namespace Tetra10 {

        constexpr int nNodes  = 10;                                        ///< Number of nodes for a Tetra10 element.
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Tetra10 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector3d& xi );
        /// @brief Evaluate Tetra10 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Eigen::Vector3i getBoundaryElementIndices( int faceID );

      } // namespace Tetra10

      namespace Hexa8 {
        constexpr int nNodes  = 8;                                         ///< Number of nodes for a Hexa8 element.
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Hexa8 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector3d& xi );
        /// @brief Evaluate Hexa8 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Eigen::Vector4i getBoundaryElementIndices( int faceID );
      } // namespace Hexa8

      namespace Hexa20 {
        constexpr int nNodes  = 20;                                        ///< Number of nodes for a Hexa20 element.
        using CoordinateSized = Eigen::Matrix< double, nNodes * nDim, 1 >; ///< Flat vector of nodal coordinates.
        using NSized          = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
        using dNdXiSized = Eigen::Matrix< double, nDim, nNodes >; ///< Matrix of shape function natural derivatives.

        /// @brief Evaluate Hexa20 shape functions at natural coordinates @p xi.
        NSized N( const Eigen::Vector3d& xi );
        /// @brief Evaluate Hexa20 shape function natural derivatives at @p xi.
        dNdXiSized dNdXi( const Eigen::Vector3d& xi );

        /// @brief Return the local node indices for boundary face @p faceID.
        Marmot::Vector8i getBoundaryElementIndices( int faceID );
      } // namespace Hexa20
    }   // namespace Spatial3D

    /**
     * @brief Boundary (surface/edge) element for computing distributed load vectors.
     *
     * Constructed from a parent element shape and a face number, this class
     * sets up the boundary quadrature points and provides methods for computing
     * scalar (pressure) and vectorial (traction) load vectors together with
     * their derivatives with respect to the parent node coordinates.
     */
    class BoundaryElement {

      /// @brief Internal data for a single quadrature point on the boundary face.
      struct BoundaryElementQuadraturePoint {
        double          weight;     ///< Quadrature weight.
        double          JxW;        ///< Jacobian determinant times quadrature weight.
        Eigen::VectorXd xi;         ///< Natural coordinates of the quadrature point.
        Eigen::VectorXd N;          ///< Shape function values at the quadrature point.
        Eigen::MatrixXd dNdXi;      ///< Shape function natural derivatives.
        Eigen::MatrixXd dx_dXi;     ///< Tangent vectors on the boundary (dx/dXi).
        Eigen::VectorXd areaVector; ///< Outward-normal area vector at the quadrature point.
      };

      const int nDim;                                                 ///< Number of spatial dimensions.

      ElementShapes boundaryShape;                                    ///< Shape of the boundary (face/edge) element.
      int           nNodes;                                           ///< Number of nodes on the boundary element.
      int           nParentCoordinates;                               ///< Total number of parent nodal coordinates.

      std::vector< BoundaryElementQuadraturePoint > quadraturePoints; ///< Quadrature points on the boundary.

      Eigen::VectorXi mapBoundaryToParentScalar;    ///< Index map: boundary scalar DOFs → parent scalar DOFs.
      Eigen::VectorXi mapBoundaryToParentVectorial; ///< Index map: boundary vectorial DOFs → parent vectorial DOFs.
      Eigen::VectorXd coordinates;                  ///< Boundary nodal coordinates (extracted from parent).

    public:
      /**
       * @brief Construct a BoundaryElement from a parent shape and face number.
       * @param[in] parentShape        Shape enum of the parent volume element.
       * @param[in] nDim               Number of spatial dimensions.
       * @param[in] parentFaceNumber   Zero-based index of the boundary face on the parent element.
       * @param[in] parentCoordinates  Flat vector of parent nodal coordinates.
       */
      BoundaryElement( ElementShapes          parentShape,
                       int                    nDim,
                       int                    parentFaceNumber,
                       const Eigen::VectorXd& parentCoordinates );

      /// @brief Compute the scalar (pressure) load vector for a unit pressure.
      Eigen::VectorXd computeScalarLoadVector();
      /// @brief Compute the derivative of the scalar load vector with respect to nodal coordinates.
      Eigen::MatrixXd computeDScalarLoadVector_dCoordinates();

      /// compute the element load vector for a unit vectorial load normal to the surface.
      Eigen::VectorXd computeSurfaceNormalVectorialLoadVector();
      /// @brief Compute the normal-traction load vector for axisymmetric elements.
      Eigen::VectorXd computeSurfaceNormalVectorialLoadVectorForAxisymmetricElements();
      /// @brief Compute the derivative of the surface-normal load vector with respect to nodal coordinates.
      Eigen::MatrixXd computeDSurfaceNormalVectorialLoadVector_dCoordinates();

      /// compute the element load vector for a unit vectorial load in a given direction.
      Eigen::VectorXd computeVectorialLoadVector( const Eigen::VectorXd& direction );
      /// @brief Compute the derivative of the directional load vector with respect to nodal coordinates.
      Eigen::MatrixXd computeDVectorialLoadVector_dCoordinates( const Eigen::VectorXd& direction );

      /// @brief Extract the boundary-node entries from a parent scalar vector.
      Eigen::VectorXd condenseParentToBoundaryScalar( const Eigen::VectorXd& parentVector );
      /// @brief Assemble a boundary scalar vector into the parent scalar vector.
      void assembleIntoParentScalar( const Eigen::VectorXd&        boundaryVector,
                                     Eigen::Ref< Eigen::VectorXd > ParentVector );
      /// @brief Assemble a boundary scalar stiffness matrix into the parent stiffness matrix.
      void assembleIntoParentStiffnessScalar( const Eigen::MatrixXd& KBoundary, Eigen::Ref< Eigen::MatrixXd > KParent );

      /// @brief Extract the boundary-node entries from a parent vectorial vector.
      Eigen::VectorXd condenseParentToBoundaryVectorial( const Eigen::VectorXd& parentVector );
      /// @brief Assemble a boundary vectorial vector into the parent vectorial vector.
      void assembleIntoParentVectorial( const Eigen::VectorXd&        boundaryVector,
                                        Eigen::Ref< Eigen::VectorXd > ParentVector );
      /// @brief Assemble a boundary vectorial stiffness matrix into the parent stiffness matrix.
      void assembleIntoParentStiffnessVectorial( const Eigen::MatrixXd&        KBoundary,
                                                 Eigen::Ref< Eigen::MatrixXd > KParent );
    };
  } // namespace FiniteElement

  namespace FiniteElement::Quadrature {
    constexpr double gp2 = 0.577350269189625764509; ///< Gauss point abscissa for 2-point rule: @f$1/\sqrt{3}@f$.
    constexpr double gp3 = 0.774596669241483;       ///< Gauss point abscissa for 3-point rule: @f$\sqrt{3/5}@f$.

    /**
     * @brief Gauss quadrature integration order options.
     */
    enum IntegrationTypes {
      FullIntegration,   ///< Full (exact) Gauss integration.
      ReducedIntegration ///< Reduced integration (one order lower).
    };

    /**
     * @brief Natural coordinates and weight for one Gauss quadrature point.
     */
    struct QuadraturePointInfo {
      Eigen::VectorXd xi;     ///< Natural coordinates of the quadrature point.
      double          weight; ///< Quadrature weight.
    };

    /**
     * @brief Return the predefined Gauss point list for a given element shape and integration type.
     * @param[in] shape           Element shape enumerator.
     * @param[in] integrationType Full or reduced integration.
     * @return Const reference to the corresponding vector of
     *         @ref Marmot::FiniteElement::Quadrature::QuadraturePointInfo entries.
     */
    const std::vector< QuadraturePointInfo >& getGaussPointInfo( Marmot::FiniteElement::ElementShapes shape,
                                                                 IntegrationTypes                     integrationType );

    /**
     * @brief Return the number of Gauss points for a given element shape and integration type.
     * @param[in] shape           Element shape enumerator.
     * @param[in] integrationType Full or reduced integration.
     * @return Number of quadrature points.
     */
    int getNumGaussPoints( Marmot::FiniteElement::ElementShapes shape, IntegrationTypes integrationType );

    namespace Spatial1D {
      constexpr int nDim = 1; ///< Number of spatial dimensions for 1-D quadrature.

      // clang-format off
            /// 1-point Gauss rule for 1-D elements.
            const std::vector< QuadraturePointInfo >  gaussPointList1 = {
                { ( Eigen::VectorXd ( 1 ) << 0 ).finished(),               2.0 }
            };

            /// 2-point Gauss rule for 1-D elements.
            const std::vector< QuadraturePointInfo >  gaussPointList2 = {
                { ( Eigen::VectorXd ( 1 ) << -gp2 ).finished(),           1.0 },
                { ( Eigen::VectorXd ( 1 ) << +gp2 ).finished(),           1.0 }
            };

            /// 3-point Gauss rule for 1-D elements.
            const std::vector< QuadraturePointInfo >  gaussPointList3 = {
                { ( Eigen::VectorXd ( 1 ) << -gp3 ).finished(),            5./9 },
                { ( Eigen::VectorXd ( 1 ) << 0.   ).finished(),            8./9 },
                { ( Eigen::VectorXd ( 1 ) << +gp3 ).finished(),            5./9 }
            };
      // clang-format on

    } // namespace Spatial1D

    namespace Spatial2D {
      constexpr int nDim = 2; ///< Number of spatial dimensions for 2-D quadrature.

      // clang-format off
            /// 1×1 Gauss rule for 2-D quadrilateral elements.
            const std::vector< QuadraturePointInfo > gaussPointList1x1 = {
                { Eigen::Vector2d::Zero(),                             4. }
            };

            /// 2×2 Gauss rule for 2-D quadrilateral elements.
            const std::vector< QuadraturePointInfo > gaussPointList2x2 = {
                { ( Eigen::Vector2d () << +gp2,     +gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << -gp2,     +gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << -gp2,     -gp2 ).finished(),   1.0 },
                { ( Eigen::Vector2d () << +gp2,     -gp2 ).finished(),   1.0 }
            };

            /// 3×3 Gauss rule for 2-D quadrilateral elements.
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

      /**
       * @brief Scale the characteristic element length in an ABAQUS-like manner.
       * @param[in,out] charElemLength Characteristic element length to be modified.
       * @param[in]     intPoint       Integration point index.
       */
      void modifyCharElemLengthAbaqusLike( double& charElemLength, int intPoint );
    } // namespace Spatial2D

    namespace Spatial3D {
      constexpr int nDim = 3; ///< Number of spatial dimensions for 3-D quadrature.

      // clang-format off
            /// 1×1×1 Gauss rule for 3-D hexahedral elements.
            const inline std::vector< QuadraturePointInfo > gaussPointList1x1x1 = {
                { Eigen::Vector3d::Zero(),                                         8.0 }
            };

            /// 1-point rule for Tetra4 elements.
            const inline std::vector< QuadraturePointInfo > gaussPointListTetra4 = {
                { (Eigen::Vector3d() << 1./4, 1./4, 1./4).finished(),  1./6}
            };

            /// 4-point rule for Tetra10 elements.
            const inline std::vector< QuadraturePointInfo > gaussPointListTetra10 = {
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20     ).finished(),  1./24},
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20,    (5+3*std::sqrt(5))/20   ).finished(),  1./24},
                { (Eigen::Vector3d() << (5-std::sqrt(5))/20,    (5+3*std::sqrt(5))/20,  (5-std::sqrt(5))/20     ).finished(),  1./24},
                { (Eigen::Vector3d() << (5+3*std::sqrt(5))/20,  (5-std::sqrt(5))/20,    (5-std::sqrt(5))/20     ).finished(),  1./24},
            };

            /// 2×2×2 Gauss rule for 3-D hexahedral elements.
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

            /// 3×3×3 Gauss rule for 3-D hexahedral elements.
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
