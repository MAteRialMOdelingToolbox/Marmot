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

using namespace Fastor;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;
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

    std::tuple< Tensor3333d, const Tensor3333d, Tensor3333d, Tensor33d > interfaceGeometrySystemCouplings(
      const Tensor33d&   N,
      const Tensor33d&   T,
      const Tensor3333d& L )
    {
      Tensor33d
        Q = Fastor::einsum< Fastor::Index< A_, i_, B_, j_ >, Fastor::Index< i_, j_ >, Fastor::OIndex< A_, B_ > >( L,
                                                                                                                  N );
      Tensor33d G = Fastor::inverse( Q );

      Tensor3333d
        A = Fastor::einsum< Fastor::Index< A_, B_ >, Fastor::Index< i_, j_ >, Fastor::OIndex< A_, i_, B_, j_ > >( G,
                                                                                                                  N );
      Tensor3333d LA  = Fastor::einsum< Fastor::Index< A_, i_, m_, n_ >,
                                       Fastor::Index< m_, n_, B_, j_ >,
                                       Fastor::OIndex< A_, i_, B_, j_ > >( L, A );
      Tensor3333d LAL = Fastor::einsum< Fastor::Index< A_, i_, m_, n_ >,
                                        Fastor::Index< m_, n_, B_, j_ >,
                                        Fastor::OIndex< A_, i_, B_, j_ > >( LA, L );
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
      Tensor3333d F                      = 1.0 * ( Fastor::einsum< Fastor::Index< A_, m_ >,
                                              Fastor::Index< m_, n_, B_, j_ >,
                                              Fastor::OIndex< A_, n_, B_, j_ > >( G_nu, L_nu ) );

      Tensor3333d Y = 1.0 * ( Fastor::einsum< Fastor::Index< A_, i_, m_, n_ >,
                                              Fastor::Index< n_, B_ >,
                                              Fastor::OIndex< A_, i_, m_, B_ > >( L_nu, G_nu ) );
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

      Tensor333d nF = Fastor::
        einsum< Fastor::Index< A_ >, Fastor::Index< i_, A_, B_, j_ >, Fastor::OIndex< i_, B_, j_ > >( normal, F );
      Tensor333d nY = Fastor::
        einsum< Fastor::Index< A_, i_, B_, j_ >, Fastor::Index< i_ >, Fastor::OIndex< A_, B_, j_ > >( Y, normal );

      Tensor333d nY_H_inv = Fastor::
        einsum< Fastor::Index< i_, j_, A_ >, Fastor::Index< A_, B_ >, Fastor::OIndex< i_, j_, B_ > >( nY, H_inv );
      Tensor333d H_inv_nF = Fastor::
        einsum< Fastor::Index< A_, B_ >, Fastor::Index< B_, i_, j_ >, Fastor::OIndex< A_, i_, j_ > >( H_inv, nF );
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< Fastor::Index< m_, i_, j_ >,
                                                Fastor::Index< m_, n_ >,
                                                Fastor::Index< n_, k_, l_ >,
                                                Fastor::OIndex< i_, j_, k_, l_ > >( nY_H_inv, G_nu, H_inv_nF );

      return std::make_tuple( B_nu, H_inv, H_inv_nF, nY_H_inv_Fn );
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Matrix6d C_nu_voigt_full = stiffnessTensor( 1.0, nu_0 );

      Tensor33d N = Fastor::einsum< Fastor::Index< i_ >, Fastor::Index< j_ >, Fastor::OIndex< i_, j_ > >( normal,
                                                                                                          normal );

      Tensor33d T = Marmot::FastorStandardTensors::Spatial3D::I - N;

      const auto  C_nu_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_nu_voigt_full );
      Tensor3333d C_nu_aibj  = convertEigenToFastor( C_nu_eigen );

      auto [Z, H_inv, H_inv_nF, nY_H_inv_Fn] = calculateMaterialMatrices( normal, N, T, C_nu_aibj );

      return { Z, H_inv, H_inv_nF, nY_H_inv_Fn };
    }

    std::tuple< Tensor3333d, Tensor33d, Tensor333d, Tensor3333d > calculateInterfaceMaterialParameters(
      const Tensor3d& normal,
      const Matrix6d& C_ep_voigt )
    {
      Tensor33d N = Fastor::einsum< Fastor::Index< i_ >, Fastor::Index< j_ >, Fastor::OIndex< i_, j_ > >( normal,
                                                                                                          normal );

      const auto  C_ep_eigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStiffness( C_ep_voigt );
      Tensor3333d C_ep_aibj  = convertEigenToFastor( C_ep_eigen );

      Tensor33d Q = Fastor::
        einsum< Fastor::Index< i_, j_, k_, l_ >, Fastor::Index< j_, l_ >, Fastor::OIndex< i_, k_ > >( C_ep_aibj, N );
      Tensor33d H_inv = Fastor::inverse( Q );

      // H_inv_nF(i,k,l) = C_ijkl * n_j  (free: i, k, l)
      Tensor333d H_inv_nF = Fastor::einsum< Fastor::Index< i_, j_, k_, l_ >,
                                            Fastor::Index< j_ >,
                                            Fastor::OIndex< i_, k_, l_ > >( C_ep_aibj, normal );
      // H_inv_Fn(i,j,k) = C_ijkl * n_l  (free: i, j, k)
      Tensor333d H_inv_Fn = Fastor::einsum< Fastor::Index< i_, j_, k_, l_ >,
                                            Fastor::Index< l_ >,
                                            Fastor::OIndex< i_, j_, k_ > >( C_ep_aibj, normal );

      // H_inv_Fn(i,j,m) * H_inv(m,r) * H_inv_nF(r,k,l) -> (i,j,k,l)
      Tensor3333d nY_H_inv_Fn = Fastor::einsum< Fastor::Index< i_, j_, m_ >,
                                                Fastor::Index< m_, o_ >,
                                                Fastor::Index< o_, k_, l_ >,
                                                Fastor::OIndex< i_, j_, k_, l_ > >( H_inv_Fn, H_inv, H_inv_nF );
      Tensor3333d Z           = C_ep_aibj - nY_H_inv_Fn;
      return { Z, Q, H_inv_nF, nY_H_inv_Fn };
    }

  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials
