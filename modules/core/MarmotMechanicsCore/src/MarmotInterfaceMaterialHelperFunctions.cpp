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
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
// #`:pragma once
// #include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/AbstractTensorFunctions.h>
#include <Fastor/tensor_algebra/einsum.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>
#include <tuple>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

enum { a, i, b, j, k, l, m, n, p, q, r, I, J };
namespace Marmot::Materials {

  namespace InterfaceMaterialHelperFunctions {
    Tensor2D compute_inv( const Tensor2D& I, Tensor2D& Q )
    {
      // Solve G in batch mode using LU decomposition
      Eigen::Matrix3d Im( I.data() );
      Eigen::Matrix3d Qm( Q.data() );

      // Compute LU decomposition with invertibility check
      Eigen::FullPivLU< Eigen::Matrix3d > lu( Qm );
      // std::cout<<"Determinant of Q: "<<lu.determinant()<<std::endl;
      if ( !lu.isInvertible() ) {
        throw std::runtime_error( "Error in compute_inv: matrix Q is singular or near-singular." );
      }

      Eigen::Matrix3d G_mat = Qm.fullPivLu().solve( Im );

      // Store result back into G (reinterpret as a 3D Fastor tensor)
      Tensor2D G( 0 );
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( G.data() ) = G_mat;

      return G;
    }
    std::tuple< Tensor4D, const Tensor4D, Tensor4D, Tensor2D > interfaceGeometrySystemCouplings( const Tensor2D& I,
                                                                                                 const Tensor2D& N,
                                                                                                 const Tensor2D& T,
                                                                                                 const Tensor4D& L )
    {
      Tensor2D Q = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i, j >, Fastor::OIndex< a, b > >( L, N );
      Tensor2D G = compute_inv( I, Q );

      Tensor4D A = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< i, j >, Fastor::OIndex< a, i, b, j > >( G, N );
      Tensor4D LA = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( L, A );
      Tensor4D LAL = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( LA, L );
      Tensor4D B = L - LAL;

      return std::make_tuple( B, L, A, G );
    }
    std::tuple< Tensor4D, Tensor4D, Tensor4D, Tensor4D, Tensor2D, Tensor4D > calculateFY( const Tensor2D& I,
                                                                                          const Tensor2D& N,
                                                                                          const Tensor2D& T,
                                                                                          const Tensor4D& C_nu_aibj )
    {
      Tensor2D G_nu;
      Tensor4D A_nu;
      Tensor4D B_nu;
      Tensor4D L_nu;

      std::tie( B_nu, L_nu, A_nu, G_nu ) = interfaceGeometrySystemCouplings( I, N, T, C_nu_aibj );
      // Tensor4D F = 1.0 * ( Fastor::einsum< Fastor::Index< a, i, m, n >,
      //                                       Fastor::Index< m, n, b, j >,
      //                                       Fastor::OIndex< a, i, b, j > >( A_nu, L_nu ) );

      // Tensor4D Y = 1.0 * ( Fastor::einsum< Fastor::Index< a, i, m, n >,
      //                                       Fastor::Index< m, n, b, j >,
      //                                       Fastor::OIndex< a, i, b, j > >( L_nu, A_nu ) );

      Tensor4D F = 1.0 * ( Fastor::einsum< Fastor::Index< a, m >,
                                           Fastor::Index< m, n, b, j >,
                                           Fastor::OIndex< a, n, b, j > >( G_nu, L_nu ) );

      Tensor4D Y = 1.0 * ( Fastor::einsum< Fastor::Index< a, i, m, n >,
                                           Fastor::Index< n, b >,
                                           Fastor::OIndex< a, i, m, b > >( L_nu, G_nu ) );
      return std::make_tuple( F, Y, A_nu, L_nu, G_nu, B_nu );
    }

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateMaterialMatrices( const Tensor1D& normal,
                                                                                    const Tensor2D& I,
                                                                                    const Tensor2D& N,
                                                                                    const Tensor2D& T,
                                                                                    const Tensor4D& C_nu_aibj )
    {

      auto [F, Y, A_nu, L_nu, G_nu, B_nu] = calculateFY( I, N, T, C_nu_aibj );

      // double H_factor = std::abs( 2.0 / E_0 - 1.0 / E_M - 1.0 / E_I );
      // double B_factor = -std::abs( E_M + E_I - 2.0 * E_0 );

      // double H_factor = std::abs( 2.0 / E_0 );
      // double B_factor = -std::abs( -2.0 * E_0 );

      // Tensor2D H = H_factor * G_nu;
      // Tensor4D Z = B_factor * B_nu;

      // std::cout<<"H_ij:\n"<<H<<std::endl;

      // std::cout<<"H_inv_ij:\n"<<std::endl;

      Tensor2D H_inv = compute_inv( I, G_nu );
      // std::cout<<"H_inv_ij:\n"<<H_inv<<std::endl;

      Tensor3D
        nF = Fastor::einsum< Fastor::Index< a >, Fastor::Index< i, a, b, j >, Fastor::OIndex< i, b, j > >( normal, F );
      Tensor3D
        nY = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i >, Fastor::OIndex< a, b, j > >( Y, normal );

      Tensor3D nY_H_inv = Fastor::
        einsum< Fastor::Index< i, j, a >, Fastor::Index< a, b >, Fastor::OIndex< i, j, b > >( nY, H_inv );
      Tensor3D
        H_inv_nF = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< b, i, j >, Fastor::OIndex< a, i, j > >( H_inv,
                                                                                                                 nF );

      // Tensor4D Yn_H_inv_Fn = Fastor::einsum< Fastor::Index< a, i, m >,
      //                                        Fastor::Index< m, n >,
      //                                        Fastor::Index< n, b, j >,
      //                                        Fastor::OIndex< a, i, b, j > >( Yn, H_inv, Fn );
      Tensor4D nY_H_inv_Fn = Fastor::einsum< Fastor::Index< m, i, j >,
                                             Fastor::Index< m, n >,
                                             Fastor::Index< n, k, l >,
                                             Fastor::OIndex< i, j, k, l > >( nY_H_inv, G_nu, H_inv_nF );

      return std::make_tuple( B_nu, H_inv, H_inv_nF, nY_H_inv_Fn );
    }

    // Convert 4th-order Fastor tensor (3x3x3x3) to Eigen 9x9 matrix
    // First two indices (ij) form rows, last two indices (kl) form columns
    Eigen::Matrix< double, 9, 9 > convert4thOrderTensorToMatrix_9x9( const Tensor4D& tensor )
    {
      Eigen::Matrix< double, 9, 9 > matrix( 9, 9 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k < 3; ++k ) {
            for ( int l = 0; l < 3; ++l ) {
              int row            = 3 * i + j; // Convert (i, j) to single index
              int col            = 3 * k + l; // Convert (k, l) to single index
              matrix( row, col ) = tensor( i, j, k, l );
            }
          }
        }
      }

      return matrix;
    }

    // Convert 3rd-order Fastor tensor (3x3x3) to Eigen 9x3 matrix
    // First two indices (ij) form rows, last index k forms columns
    Eigen::Matrix< double, 9, 3 > convert3rdOrderTensorToMatrix_9x3( const Tensor3D& tensor )
    {
      Eigen::Matrix< double, 9, 3 > matrix( 9, 3 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k < 3; ++k ) {
            int row          = 3 * i + j; // Convert (i, j) to single index
            matrix( row, k ) = tensor( i, j, k );
          }
        }
      }

      return matrix;
    }

    // Convert 3rd-order Fastor tensor (3x3x3) to Eigen 3x9 matrix
    // First index i forms rows, last two indices (jk) form columns
    Eigen::Matrix< double, 3, 9 > convert3rdOrderTensorToMatrix_3x9( const Tensor3D& tensor )
    {
      Eigen::Matrix< double, 3, 9 > matrix( 3, 9 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k < 3; ++k ) {
            int col          = 3 * j + k; // Convert (j, k) to single index
            matrix( i, col ) = tensor( i, j, k );
          }
        }
      }

      return matrix;
    }

    // Convert 2nd-order Fastor tensor (3x3) to Eigen 3x3 matrix
    Eigen::Matrix< double, 3, 3 > convert2ndOrderTensorToMatrix_3x3( const Tensor2D& tensor )
    {
      Eigen::Matrix< double, 3, 3 > matrix( 3, 3 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          matrix( i, j ) = tensor( i, j );
        }
      }
      return matrix;
    }

    Tensor4D voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness )
    {
      using namespace Marmot::ContinuumMechanics::TensorUtility::IndexNotation;

      Fastor::Tensor< double, 3, 3, 3, 3 > stiffness;
      stiffness.zeros(); // Set to zero

      int row, col;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j ) {
          row = toVoigt< 3 >( i, j );
          for ( int k = 0; k < 3; ++k )
            for ( int l = 0; l < 3; ++l ) {
              col                     = toVoigt< 3 >( k, l );
              stiffness( i, j, k, l ) = voigtStiffness( row, col );
            }
        }
      return stiffness;
    }

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateInterfaceMaterialParameters( const Tensor1D& normal,
                                                                                               const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Eigen::Matrix< double, 6, 6 > C_nu_voigt_full = stiffnessTensor( 1.0, nu_0 );

      Tensor2D I = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

      Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

      Tensor2D T = I - N;

      Tensor4D C_nu_aibj = voigtToStiffness( C_nu_voigt_full );

      auto [Z, H_inv, H_inv_nF, nY_H_inv_Fn] = calculateMaterialMatrices( normal, I, N, T, C_nu_aibj );

      return { Z, H_inv, H_inv_nF, nY_H_inv_Fn };
    }

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateInterfaceMaterialParameters(
      const Tensor1D&                      normal,
      const Eigen::Matrix< double, 6, 6 >& C_ep_voigt )
    {
      Tensor2D I = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

      Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );


      Tensor4D C_ep_aibj = voigtToStiffness( C_ep_voigt );

      Tensor2D
        Q = Fastor::einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j, l >, Fastor::OIndex< i, k > >( C_ep_aibj,
                                                                                                          N );
      Tensor2D H_inv = compute_inv( I, Q );

      // H_inv_nF(i,k,l) = C_ijkl * n_j  (free: i, k, l)
      Tensor3D H_inv_nF = Fastor::
        einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j >, Fastor::OIndex< i, k, l > >( C_ep_aibj, normal );
      // H_inv_Fn(i,j,k) = C_ijkl * n_l  (free: i, j, k)
      Tensor3D H_inv_Fn = Fastor::
        einsum< Fastor::Index< i, j, k, l >, Fastor::Index< l >, Fastor::OIndex< i, j, k > >( C_ep_aibj, normal );

      // H_inv_Fn(i,j,m) * H_inv(m,r) * H_inv_nF(r,k,l) -> (i,j,k,l)
      Tensor4D nY_H_inv_Fn = Fastor::einsum< Fastor::Index< i, j, m >,
                                             Fastor::Index< m, r >,
                                             Fastor::Index< r, k, l >,
                                             Fastor::OIndex< i, j, k, l > >( H_inv_Fn, H_inv, H_inv_nF );
      Tensor4D Z           = C_ep_aibj - nY_H_inv_Fn;
      return { Z, Q, H_inv_nF, nY_H_inv_Fn };
    }

  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials