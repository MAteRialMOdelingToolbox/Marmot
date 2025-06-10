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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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

#include "Marmot/MarmotFastorTensorBasics.h"
#include <Fastor/tensor/Tensor.h>

namespace Marmot::ContinuumMechanics {

  namespace StressMeasures {

    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    template < typename T >
    Tensor33t< T > KirchhoffStressFromPK2( const Tensor33t< T >& PK2, const Tensor33t< T >& F )
    {
      // tau = F * PK2 * F^T  = F_iI * S_IJ * F_Jj
      const Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

      return tau;
    }

    namespace FirstOrderDerived {

      template < typename T >
      std::tuple< Tensor33t< T >, Tensor3333t< T >, Tensor3333t< T > > KirchhoffStressFromPK2(
        const Tensor33t< T >& PK2,
        const Tensor33t< T >& F )
      {
        const auto&          I   = Spatial3D::I;
        const Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );
        /* const Tensor3333t< T > dTau_dPK2 = einsum< iK, jL, to_ijKL >( einsum< iI, IK >( F, I ), */
        /*                                                               transpose( einsum< JL, jJ >( I, F ) ) ); */
        const Tensor3333t< T > dTau_dPK2 = einsum< iI, jJ, to_ijIJ >( F, F );
        /* const Tensor3333t< T > dTau_dF   = einsum< iK, IL, IJ, jJ, to_ijKL >( I, I, PK2, F ) + */
        /*                                  einsum< iI, IJ, jK, JL, to_ijKL >( F, PK2, I, I ); */

        const Tensor33t< T >   S_F     = einsum< KJ, jJ >( PK2, F );
        const Tensor3333t< T > dTau_dF = einsum< ik, jK, to_ijkK >( I, transpose( S_F ) ) +
                                         einsum< Ki, jk, to_ijkK >( S_F, I );

        return { tau, dTau_dPK2, dTau_dF };
      }
    } // namespace FirstOrderDerived

  }   // namespace StressMeasures

} // namespace Marmot::ContinuumMechanics
