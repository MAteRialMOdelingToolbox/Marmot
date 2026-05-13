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

#include "Marmot/MarmotFastorTensorBasics.h"

using namespace Eigen;
using namespace Fastor;
using Marmot::FastorStandardTensors::Tensor3333d;
using Marmot::FastorStandardTensors::Tensor333d;
using Marmot::FastorStandardTensors::Tensor33d;
using Marmot::FastorStandardTensors::Tensor3d;

enum { a, i, b, j, k, l, m, n, p, q, r, I, J };
namespace Marmot::Materials {
  namespace InterfaceMaterialHelperFunctions {

    Tensor3333d convertEigenToFastor( const Marmot::EigenTensors::Tensor3333d& tensorEigen )
    {
      Tensor3333d tensorFastor;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          for ( int k = 0; k < 3; ++k )
            for ( int l = 0; l < 3; ++l )
              tensorFastor( i, j, k, l ) = tensorEigen( i, j, k, l );

      return tensorFastor;
    }

    Tensor33d compute_inv( const Tensor33d& Q )
    {
      return Fastor::inverse( Q );
    }
    std::tuple< Tensor3333d, const Tensor3333d, Tensor3333d, Tensor33d > interfaceGeometrySystemCouplings(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& L )
    {
      Tensor33d Q = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i, j >, Fastor::OIndex< a, b > >( L,
                                                                                                                  N );
      Tensor33d G = compute_inv( Q );

      Tensor3333d A  = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< i, j >, Fastor::OIndex< a, i, b, j > >( G,
                                                                                                                    N );
      Tensor3333d LA = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( L, A );
      Tensor3333d LAL = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( LA, L );
      Tensor3333d B = L - LAL;

      return std::make_tuple( B, L, A, G );
    }
    std::tuple< Tensor3333d, Tensor3333d, Tensor3333d, Tensor3333d, Tensor33d, Tensor3333d > calculateFY(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& C_nu_aibj )
    {
      Tensor33d   G_nu;
      Tensor3333d A_nu;
      Tensor3333d B_nu;
      Tensor3333d L_nu;

      std::tie( B_nu, L_nu, A_nu, G_nu ) = interfaceGeometrySystemCouplings( N, T, C_nu_aibj );
      Tensor3333d F                      = 1.0 * ( Fastor::einsum< Fastor::Index< a, m >,
                                              Fastor::Index< m, n, b, j >,
                                              Fastor::OIndex< a, n, b, j > >( G_nu, L_nu ) );

      Tensor3333d Y = 1.0 * ( Fastor::einsum< Fastor::Index< a, i, m, n >,
                                              Fastor::Index< n, b >,
                                              Fastor::OIndex< a, i, m, b > >( L_nu, G_nu ) );
      return std::make_tuple( F, Y, A_nu, L_nu, G_nu, B_nu );
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateMaterialMatrices(
      const Tensor3d&    normal,
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& C_nu_aibj )
    {

      auto [F, Y, A_nu, L_nu, G_nu, B_nu] = calculateFY( N, T, C_nu_aibj );

      Tensor33d H_inv = compute_inv( G_nu );

      Tensor333d
        nF = Fastor::einsum< Fastor::Index< a >, Fastor::Index< i, a, b, j >, Fastor::OIndex< i, b, j > >( normal, F );
      Tensor333d
        nY = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i >, Fastor::OIndex< a, b, j > >( Y, normal );

      Tensor333d nY_H_inv = Fastor::
        einsum< Fastor::Index< i, j, a >, Fastor::Index< a, b >, Fastor::OIndex< i, j, b > >( nY, H_inv );
      Tensor333d
        H_inv_nF = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< b, i, j >, Fastor::OIndex< a, i, j > >( H_inv,
                                                                                                                 nF );
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< Fastor::Index< m, i, j >,
                                                Fastor::Index< m, n >,
                                                Fastor::Index< n, k, l >,
                                                Fastor::OIndex< i, j, k, l > >( nY_H_inv, G_nu, H_inv_nF );

      return std::make_tuple( B_nu, H_inv, H_inv_nF, nY_H_inv_Fn );
    }

    Eigen::Matrix< double, 9, 9 > convert4thOrderTensorToMatrix_9x9( const Tensor3333d& tensor )
    {
      return Marmot::ContinuumMechanics::TensorUtility::convert4thOrderTensorToMatrix_9x9( tensor );
    }

    Eigen::Matrix< double, 9, 3 > convert3rdOrderTensorToMatrix_9x3( const Tensor333d& tensor )
    {
      return Marmot::ContinuumMechanics::TensorUtility::convert3rdOrderTensorToMatrix_9x3( tensor );
    }

    Eigen::Matrix< double, 3, 9 > convert3rdOrderTensorToMatrix_3x9( const Tensor333d& tensor )
    {
      return Marmot::ContinuumMechanics::TensorUtility::convert3rdOrderTensorToMatrix_3x9( tensor );
    }

    Eigen::Matrix< double, 3, 3 > convert2ndOrderTensorToMatrix_3x3( const Tensor33d& tensor )
    {
      return Marmot::ContinuumMechanics::TensorUtility::convert2ndOrderTensorToMatrix_3x3( tensor );
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Eigen::Matrix< double, 6, 6 > C_nu_voigt_full = stiffnessTensor( 1.0, nu_0 );

      Tensor33d N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

      Tensor33d T = Marmot::FastorStandardTensors::Spatial3D::I - N;

      const auto  C_nu_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_nu_voigt_full );
      Tensor3333d C_nu_aibj  = convertEigenToFastor( C_nu_eigen );

      auto [Z, H_inv, H_inv_nF, nY_H_inv_Fn] = calculateMaterialMatrices( normal, N, T, C_nu_aibj );

      return { Z, H_inv, H_inv_nF, nY_H_inv_Fn };
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d&                      normal,
      const Eigen::Matrix< double, 6, 6 >& C_ep_voigt )
    {
      Tensor33d N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

      const auto  C_ep_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_ep_voigt );
      Tensor3333d C_ep_aibj  = convertEigenToFastor( C_ep_eigen );

      Tensor33d
        Q = Fastor::einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j, l >, Fastor::OIndex< i, k > >( C_ep_aibj,
                                                                                                          N );
      Tensor33d H_inv = compute_inv( Q );

      // H_inv_nF(i,k,l) = C_ijkl * n_j  (free: i, k, l)
      Tensor333d H_inv_nF = Fastor::
        einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j >, Fastor::OIndex< i, k, l > >( C_ep_aibj, normal );
      // H_inv_Fn(i,j,k) = C_ijkl * n_l  (free: i, j, k)
      Tensor333d H_inv_Fn = Fastor::
        einsum< Fastor::Index< i, j, k, l >, Fastor::Index< l >, Fastor::OIndex< i, j, k > >( C_ep_aibj, normal );

      // H_inv_Fn(i,j,m) * H_inv(m,r) * H_inv_nF(r,k,l) -> (i,j,k,l)
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< Fastor::Index< i, j, m >,
                                                Fastor::Index< m, r >,
                                                Fastor::Index< r, k, l >,
                                                Fastor::OIndex< i, j, k, l > >( H_inv_Fn, H_inv, H_inv_nF );
      Tensor3333d Z           = C_ep_aibj - nY_H_inv_Fn;
      return { Z, Q, H_inv_nF, nY_H_inv_Fn };
    }

  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
