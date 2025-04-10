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

namespace Marmot::ContinuumMechanics {

  namespace DeformationMeasures {

    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    template < typename T >
    Tensor33t< T > CauchyGreen( const Tensor33t< T > F_ )
    {
      // C = F ^ T * F; C_IJ = F_iI F_iJ
      Tensor33t< T > C = einsum< iI, iJ >( F_, F_ );
      return C;
    }

    namespace FirstOrderDerived {

      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > CauchyGreen( const Tensor33t< T > F )
      {
        Tensor33t< T > C = einsum< iI, iJ >( F, F );

        // dC_IJ / dF_F_KL = dF_iI / dF_KL * F_iJ + F_iI * dF_iJ / dF_F_KL
        const Tensor3333t< T > dC_dF = einsum< IK, kJ, to_IJkK >( Spatial3D::I, F ) +
                                       einsum< kI, JK, to_IJkK >( F, Spatial3D::I );

        return { C, dC_dF };
      }

    } // namespace FirstOrderDerived

  }   // namespace DeformationMeasures
} // namespace Marmot::ContinuumMechanics
