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

#include "Marmot/MarmotMicromorphicTensorBasics.h"

namespace Marmot::ContinuumMechanics {

  namespace DeformationMeasures {

    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    template < typename T >
    Tensor33t< T > CauchyGreen( const Tensor33t< T > F_ )
    {
      Tensor33t< T > C = einsum< KI, KJ >( F_, F_ );
      return C;
    }

    namespace FirstOrderDerived {

      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > CauchyGreen( const Tensor33t< T > F_ )
      {
        Tensor33t< T > C = einsum< KI, KJ >( F_, F_ );

        const auto&            I     = FastorStandardTensors::Spatial3D::I;
        const Tensor3333t< T > dC_dF = einsum< IL, KJ, to_IJKL >( I, F_ ) + einsum< JL, KI, to_IJKL >( I, F_ );

        return { C, dC_dF };
      }

    } // namespace FirstOrderDerived

  } // namespace DeformationMeasures
} // namespace Marmot::ContinuumMechanics
