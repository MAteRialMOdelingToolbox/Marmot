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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensorExponential.h"
#include <cmath>

/**
 * @file MarmotFiniteStrainPlasticity.h
 * @brief Helpers for integrating the plastic flow rule in finite-strain plasticity.
 *
 * Provides exponential-map and first-order approximations for the incremental
 * plastic deformation gradient together with their derivatives with respect to
 * the plastic velocity gradient.
 */

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

      /** Computes the incremental plastic deformation gradient from the plastic velocity gradient
       *  using the exponential map, evaluated by scaling and squaring.
       *
       *  The result is that of exponentialMap(), but the series is only evaluated for
       *  \f$ \Delta\boldsymbol{G}^p / 2^s \f$ with \f$ |\Delta\boldsymbol{G}^p| / 2^s \leq 1/2 \f$, and then squared
       *  \f$ s \f$ times,
       *  \f[
       *    \exp( \Delta\boldsymbol{G}^p ) = \left( \exp( \Delta\boldsymbol{G}^p / 2^s ) \right)^{2^s}.
       *  \f]
       *  For \f$ |\Delta\boldsymbol{G}^p| \leq 1/2 \f$, \f$ s = 0 \f$ and the result is identical to that of
       *  exponentialMap(). Larger increments, e.g. those of the intermediate iterates of a return mapping, are
       *  represented as well, where the truncated series of exponentialMap() fails.
       *
       *  @tparam T Scalar type, real or complex (the scaling is chosen from the real part).
       *  @param dGp Incremental plastic velocity gradient.
       *  @return Incremental plastic deformation gradient.
       */
      template < typename T >
      Tensor33t< T > exponentialMapScalingAndSquaring( const Tensor33t< T >& dGp )
      {
        double norm2 = 0.0;
        for ( int i = 0; i < 9; i++ )
          norm2 += std::pow( std::abs( Math::makeReal( dGp.data()[i] ) ), 2 );
        const int      s   = norm2 > 0.25 ? int( std::ceil( std::log2( std::sqrt( norm2 ) / 0.5 ) ) ) : 0;
        Tensor33t< T > dFp = exponentialMap(
          Tensor33t< T >( multiplyFastorTensorWithScalar( dGp, T( std::ldexp( 1.0, -s ) ) ) ) );
        for ( int k = 0; k < s; k++ )
          dFp = Tensor33t< T >( dFp % dFp );
        return dFp;
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
