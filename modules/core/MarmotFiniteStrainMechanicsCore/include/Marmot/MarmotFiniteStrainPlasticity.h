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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include "Marmot/MarmotTensorExponential.h"

namespace Marmot {
  namespace ContinuumMechanics::FiniteStrain::Plasticity {

    namespace FlowIntegration {

      using namespace FastorStandardTensors;
      using namespace Fastor;
      template < typename T >
      Tensor33t< T > exponentialMap( const Tensor33t< T >& dGp )
      {
        const Tensor33t< T > dFpT = TensorUtility::TensorExponential::computeTensorExponential( dGp, 15, 1e-14, 1e-14 );
        const Tensor33t< T > out  = permute< Index< 1, 0 > >( dFpT );
        return out;
      }
      namespace FirstOrderDerived {
        std::pair< Tensor33d, Tensor3333d > explicitIntegration( const Tensor33d& deltaGp );

        std::pair< Tensor33d, Tensor3333d > exponentialMap( const Tensor33d& deltaGp );
      } // namespace FirstOrderDerived

    }   // namespace FlowIntegration
  }     // namespace ContinuumMechanics::FiniteStrain::Plasticity
} // namespace Marmot
