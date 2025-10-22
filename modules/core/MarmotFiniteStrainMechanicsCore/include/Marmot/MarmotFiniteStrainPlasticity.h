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

      /** Computes the incremental plastic deformation gradient from the plastic velocity gradient
       *  using the exponential map.
       *
       *  The exponential map is given by
       *  \f[
       *    \Delta F_{\bar I J} = \exp\left(\Delta \lambda \frac{\partial g}{\partial M_{\bar K L}}\right)_{J \bar{I}}
       *  \f]
       *  where \f$ \Delta \lambda \frac{\partial g}{\partial M_{\bar K L}} \f$ is the incremental plastic velocity
       * gradient.
       *
       *  @tparam T Scalar type.
       *  @param dGp Incremental plastic velocity gradient.
       *  @return Incremental plastic deformation gradient.
       */
      template < typename T >
      Tensor33t< T > exponentialMap( const Tensor33t< T >& dGp )
      {
        const Tensor33t< T > dFpT = TensorUtility::TensorExponential::computeTensorExponential( dGp, 15, 1e-14, 1e-14 );
        const Tensor33t< T > out  = permute< Index< 1, 0 > >( dFpT );
        return out;
      }
      namespace FirstOrderDerived {

        /** Computes the incremental plastic deformation gradient from the plastic velocity gradient
         *  using a first order approximation.
         *  The first order approximation is given by
         *  \f[
         *    \Delta F_{\bar I J} = \left( \delta_{\bar K L} + \Delta \lambda \frac{\partial f}{\partial M_{\bar K
         * L}}\right)_{J \bar{I}} \f] where \f$ \Delta \lambda \frac{\partial f}{\partial M_{\bar K L}} \f$ is the
         * incremental plastic velocity gradient.
         *
         *  @param deltaGp Incremental plastic velocity gradient.
         *  @return A pair of the incremental plastic deformation gradient and its derivative w.r.t. the plastic
         * velocity gradient.
         */
        std::pair< Tensor33d, Tensor3333d > explicitIntegration( const Tensor33d& deltaGp );

        /** Computes the incremental plastic deformation gradient from the plastic velocity gradient
         *  using the exponential map.
         *  The exponential map is given by
         *  \f[
         *    \Delta F_{\bar I J} = \exp\left(\Delta \lambda \frac{\partial f}{\partial M_{\bar I J}}\right)_{J \bar{I}}
         *  \f]
         *  where \f$ \Delta \lambda \frac{\partial f}{\partial M_{\bar I J}} \f$ is the incremental plastic velocity
         * gradient.
         *
         *  @param deltaGp Incremental plastic velocity gradient.
         *  @return A pair of the incremental plastic deformation gradient and its derivative w.r.t. the plastic
         * velocity gradient.
         */
        std::pair< Tensor33d, Tensor3333d > exponentialMap( const Tensor33d& deltaGp );
      } // namespace FirstOrderDerived

    }   // namespace FlowIntegration
  }     // namespace ContinuumMechanics::FiniteStrain::Plasticity
} // namespace Marmot
