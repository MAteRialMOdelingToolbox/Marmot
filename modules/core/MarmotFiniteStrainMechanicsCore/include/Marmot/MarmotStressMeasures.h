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

    /**
     * @brief Computes the Kirchhoff stress from the 2nd Piola-Kirchhoff stress and the deformation gradient.
     *
     * The Kirchhoff stress is computed as
     * \f[
     *   \boldsymbol{\tau} = \boldsymbol{F}  \boldsymbol{S}  \boldsymbol{F}^T
     * \f]
     * or in index notation
     * \f[
     * \tau_{ij} = F_{iI} S_{IJ} F_{jJ}.
     * \f]
     *
     * @tparam T Scalar type, e.g. double, float
     * @param PK2 2nd Piola-Kirchhoff stress
     * @param F Deformation gradient
     * @return Kirchhoff stress
     */
    template < typename T >
    Tensor33t< T > KirchhoffStressFromPK2( const Tensor33t< T >& PK2, const Tensor33t< T >& F )
    {
      const Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

      return tau;
    }

    namespace FirstOrderDerived {

      /** @brief Computes the Kirchhoff stress from the 2nd Piola-Kirchhoff stress and the deformation gradient and the
       * respective partial derivatives.
       *
       * The Kirchoff stress is computed as
       * \f[
       *   \boldsymbol{\tau} = \boldsymbol{F}  \boldsymbol{S}  \boldsymbol{F}^T
       * \f]
       * or in index notation
       * \f[
       * \tau_{ij} = F_{iI} S_{IJ} F_{jJ}.
       * \f]
       * Additionally, the derivatives with respect to \f$ \boldsymbol{S} \f$ and \f$ \boldsymbol{F} \f$ are computed
       * in index notation as \f[ \frac{\partial \tau_{ij}}{\partial S_{IJ}} = F_{iI} F_{jJ} \f] and \f[ \frac{\partial
       * \tau_{ij}}{\partial F_{kL}} = \delta_{ik} S_{LJ} F_{jJ} + F_{iI} S_{IL} \delta_{jk} \f]
       *
       * @tparam T Scalar type, e.g. double, float
       * @param PK2 2nd Piola-Kirchhoff stress
       * @param F Deformation gradient
       * @return A tuple containing Kirchhoff stress and its derivatives w.r.t \f$\boldsymbol{S}\f$ and
       * \f$\boldsymbol{F}\f$
       */
      template < typename T >
      std::tuple< Tensor33t< T >, Tensor3333t< T >, Tensor3333t< T > > KirchhoffStressFromPK2(
        const Tensor33t< T >& PK2,
        const Tensor33t< T >& F )
      {
        const auto&          I   = Spatial3D::I;
        const Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

        const Tensor3333t< T > dTau_dPK2 = einsum< iI, jJ, to_ijIJ >( F, F );

        const Tensor33t< T >   S_F     = einsum< KJ, jJ >( PK2, F );
        const Tensor3333t< T > dTau_dF = einsum< ik, jK, to_ijkK >( I, transpose( S_F ) ) +
                                         einsum< Ki, jk, to_ijkK >( S_F, I );

        return { tau, dTau_dPK2, dTau_dF };
      }
    } // namespace FirstOrderDerived

  }   // namespace StressMeasures

} // namespace Marmot::ContinuumMechanics
