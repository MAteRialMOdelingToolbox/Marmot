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
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotElasticity.h"

using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;
namespace Marmot::Materials {
  namespace InterfaceMaterialHelperFunctions {

    std::tuple< Tensor3333d, const Tensor3333d, Tensor3333d, Tensor33d > interfaceGeometrySystemCouplings(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& L )
    {
      Tensor33d Q = Fastor::einsum< AiBj, ij, to_AB >( L, N );
      Tensor33d G = Fastor::inverse( Q );

      Tensor3333d A   = Fastor::einsum< AB, ij, to_AiBj >( G, N );
      Tensor3333d LA  = Fastor::einsum< Aimn, mnBj, to_AiBj >( L, A );
      Tensor3333d LAL = Fastor::einsum< Aimn, mnBj, to_AiBj >( LA, L );
      Tensor3333d B   = L - LAL;

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
      Tensor3333d F                      = Fastor::einsum< Am, mnBj, to_AnBj >( G_nu, L_nu );
      Tensor3333d Y                      = Fastor::einsum< Aimn, nB, to_AimB >( L_nu, G_nu );
      return std::make_tuple( F, Y, A_nu, L_nu, G_nu, B_nu );
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateMaterialMatrices(
      const Tensor3d&    normal,
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& C_nu_aibj )
    {

      auto [F, Y, A_nu, L_nu, G_nu, B_nu] = calculateFY( N, T, C_nu_aibj );

      Tensor33d H_inv = Fastor::inverse( G_nu );

      Tensor333d  nF          = Fastor::einsum< A, iABj, to_iBj >( normal, F );
      Tensor333d  nY          = Fastor::einsum< AiBj, i, to_ABj >( Y, normal );
      Tensor333d  nY_H_inv    = Fastor::einsum< ijA, AB, to_ijB >( nY, H_inv );
      Tensor333d  H_inv_nF    = Fastor::einsum< AB, Bij, to_Aij >( H_inv, nF );
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< mij, mn, nkl, to_ijkl >( nY_H_inv, G_nu, H_inv_nF );

      return std::make_tuple( B_nu, H_inv, H_inv_nF, nY_H_inv_Fn );
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Matrix6d C_nu_voigt_full = stiffnessTensor( 1.0, nu_0 );

      Tensor33d N = Fastor::einsum< i, j, to_ij >( normal, normal );

      Tensor33d T = Marmot::FastorStandardTensors::Spatial3D::I - N;

      const auto  C_nu_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_nu_voigt_full );
      Tensor3333d C_nu_aibj( C_nu_eigen.data(), Fastor::ColumnMajor );

      auto [Z, H_inv, H_inv_nF, nY_H_inv_Fn] = calculateMaterialMatrices( normal, N, T, C_nu_aibj );

      return { Z, H_inv, H_inv_nF, nY_H_inv_Fn };
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const Matrix6d& C_ep_voigt )
    {
      Tensor33d N = Fastor::einsum< i, j, to_ij >( normal, normal );

      const auto  C_ep_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_ep_voigt );
      Tensor3333d C_ep_aibj( C_ep_eigen.data(), Fastor::ColumnMajor );

      Tensor33d Q     = Fastor::einsum< ijkl, jl, to_ik >( C_ep_aibj, N );
      Tensor33d H_inv = Fastor::inverse( Q );

      // H_inv_nF(i,k,l) = C_ijkl * n_j  (free: i, k, l)
      Tensor333d H_inv_nF = Fastor::einsum< ijkl, j, to_ikl >( C_ep_aibj, normal );
      // H_inv_Fn(i,j,k) = C_ijkl * n_l  (free: i, j, k)
      Tensor333d H_inv_Fn = Fastor::einsum< ijkl, l, to_ijk >( C_ep_aibj, normal );

      // H_inv_Fn(i,j,m) * H_inv(m,r) * H_inv_nF(r,k,l) -> (i,j,k,l)
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< ijm, mo, okl, to_ijkl >( H_inv_Fn, H_inv, H_inv_nF );
      Tensor3333d Z           = C_ep_aibj - nY_H_inv_Fn;
      return { Z, Q, H_inv_nF, nY_H_inv_Fn };
    }

  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
