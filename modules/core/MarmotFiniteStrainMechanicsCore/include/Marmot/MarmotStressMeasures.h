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
#include <Fastor/tensor/Tensor.h>
#include <tuple>
#include <utility>

namespace Marmot::ContinuumMechanics::StressMeasures {

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
  Math::Tensors::FastorTensors::Tensor33t< T > KirchhoffStressFromPK2(
    const Math::Tensors::FastorTensors::Tensor33t< T >& PK2,
    const Math::Tensors::FastorTensors::Tensor33t< T >& F )
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;
    const Math::Tensors::FastorTensors::Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

    return tau;
  }

  namespace FirstOrder {

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
    std::tuple< Math::Tensors::FastorTensors::Tensor33t< T >,
                Math::Tensors::FastorTensors::Tensor3333t< T >,
                Math::Tensors::FastorTensors::Tensor3333t< T > >
    KirchhoffStressFromPK2( const Math::Tensors::FastorTensors::Tensor33t< T >& PK2,
                            const Math::Tensors::FastorTensors::Tensor33t< T >& F )
    {
      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;
      const auto&                                        I   = Math::Tensors::FastorTensors::Spatial3D::I;
      const Math::Tensors::FastorTensors::Tensor33t< T > tau = einsum< iI, IJ, jJ, to_ij >( F, PK2, F );

      const Math::Tensors::FastorTensors::Tensor3333t< T > dTau_dPK2 = einsum< iI, jJ, to_ijIJ >( F, F );

      const Math::Tensors::FastorTensors::Tensor33t< T >   S_F     = einsum< KJ, jJ >( PK2, F );
      const Math::Tensors::FastorTensors::Tensor3333t< T > dTau_dF = einsum< ik, jK, to_ijkK >( I, transpose( S_F ) ) +
                                                                     einsum< Ki, jk, to_ijkK >( S_F, I );

      return { tau, dTau_dPK2, dTau_dF };
    }
    /** @brief Computes the Kirchhoff stress from the 2nd Piola-Kirchhoff stress and its total derivative with respect
     * to the deformation gradient, for a 2nd Piola-Kirchhoff stress that depends on \f$\boldsymbol{F}\f$ through the
     * right Cauchy-Green tensor \f$\boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F}\f$.
     *
     * With \f$\tau_{ij} = F_{iI} S_{IJ} F_{jJ}\f$ and \f$\boldsymbol{S} = \boldsymbol{S}(\boldsymbol{C})\f$,
     * \f[
     *   \frac{d \tau_{ij}}{d F_{kM}} = \delta_{ik} S_{MJ} F_{jJ} + F_{iI} S_{IM} \delta_{jk}
     *   + F_{iI} F_{jJ} \left( \frac{\partial S_{IJ}}{\partial C_{ML}} + \frac{\partial S_{IJ}}{\partial C_{LM}}
     *   \right) F_{kL}.
     * \f]
     * This is the tangent
     * \f$\partial\boldsymbol\tau/\partial\boldsymbol S : \partial\boldsymbol S/\partial\boldsymbol C :
     * \partial\boldsymbol C/\partial\boldsymbol F + \partial\boldsymbol\tau/\partial\boldsymbol F\f$ assembled
     * from #KirchhoffStressFromPK2 and DeformationMeasures::FirstOrder::rightCauchyGreen, but evaluated without the
     * fourth-order tensors \f$\partial\boldsymbol\tau/\partial\boldsymbol S\f$ and
     * \f$\partial\boldsymbol C/\partial\boldsymbol F\f$ and their two double contractions: the deformation
     * gradient is applied index by index. Neither \f$\boldsymbol S\f$ nor
     * \f$\partial\boldsymbol S/\partial\boldsymbol C\f$ needs to be symmetric.
     *
     * For a hyperelastic potential \f$\Psi(\boldsymbol C)\f$, pass \f$\boldsymbol S = 2\,\partial\Psi/\partial
     * \boldsymbol C\f$ and \f$\partial\boldsymbol S/\partial\boldsymbol C = 2\,\partial^2\Psi/\partial\boldsymbol
     * C\partial\boldsymbol C\f$; for a derivative w.r.t. the Green-Lagrange strain, \f$\partial\boldsymbol
     * S/\partial\boldsymbol C = \tfrac{1}{2}\,\partial\boldsymbol S/\partial\boldsymbol E\f$.
     *
     * @tparam T Scalar type, e.g. double, float
     * @param PK2 2nd Piola-Kirchhoff stress
     * @param dPK2_dC Derivative of the 2nd Piola-Kirchhoff stress w.r.t. the right Cauchy-Green tensor
     * @param F Deformation gradient
     * @return A pair containing the Kirchhoff stress and its total derivative w.r.t. \f$\boldsymbol{F}\f$
     */
    template < typename T >
    std::pair< Math::Tensors::FastorTensors::Tensor33t< T >, Math::Tensors::FastorTensors::Tensor3333t< T > > KirchhoffStressAndTangentFromPK2(
      const Math::Tensors::FastorTensors::Tensor33t< T >&   PK2,
      const Math::Tensors::FastorTensors::Tensor3333t< T >& dPK2_dC,
      const Math::Tensors::FastorTensors::Tensor33t< T >&   F )
    {
      using Math::Tensors::FastorTensors::Tensor3333t;
      using Math::Tensors::FastorTensors::Tensor33t;

      std::pair< Tensor33t< T >, Tensor3333t< T > > result;
      auto& [tau, dTau_dF] = result;

      // FS_iJ = F_iI S_IJ and SF_Kj = S_KJ F_jJ
      Tensor33t< T > FS, SF;
      for ( int i = 0; i < 3; i++ )
        for ( int J = 0; J < 3; J++ ) {
          FS( i, J ) = F( i, 0 ) * PK2( 0, J ) + F( i, 1 ) * PK2( 1, J ) + F( i, 2 ) * PK2( 2, J );
          SF( J, i ) = PK2( J, 0 ) * F( i, 0 ) + PK2( J, 1 ) * F( i, 1 ) + PK2( J, 2 ) * F( i, 2 );
        }

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          tau( i, j ) = FS( i, 0 ) * F( j, 0 ) + FS( i, 1 ) * F( j, 1 ) + FS( i, 2 ) * F( j, 2 );

      // A_ijKL = F_iI F_jJ dS_IJ/dC_KL, one index at a time
      Tensor3333t< T > B, A;
      for ( int i = 0; i < 3; i++ )
        for ( int J = 0; J < 3; J++ )
          for ( int K = 0; K < 3; K++ )
            for ( int L = 0; L < 3; L++ )
              B( i, J, K, L ) = F( i, 0 ) * dPK2_dC( 0, J, K, L ) + F( i, 1 ) * dPK2_dC( 1, J, K, L ) +
                                F( i, 2 ) * dPK2_dC( 2, J, K, L );

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int K = 0; K < 3; K++ )
            for ( int L = 0; L < 3; L++ )
              A( i, j, K, L ) = F( j, 0 ) * B( i, 0, K, L ) + F( j, 1 ) * B( i, 1, K, L ) + F( j, 2 ) * B( i, 2, K, L );

      // dC_KL/dF_kM = delta_KM F_kL + F_kK delta_LM
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          for ( int k = 0; k < 3; k++ )
            for ( int M = 0; M < 3; M++ ) {
              T dTau = ( A( i, j, M, 0 ) + A( i, j, 0, M ) ) * F( k, 0 ) +
                       ( A( i, j, M, 1 ) + A( i, j, 1, M ) ) * F( k, 1 ) +
                       ( A( i, j, M, 2 ) + A( i, j, 2, M ) ) * F( k, 2 );
              if ( i == k )
                dTau += SF( M, j );
              if ( j == k )
                dTau += FS( i, M );
              dTau_dF( i, j, k, M ) = dTau;
            }

      return result;
    }
  } // namespace FirstOrder

} // namespace Marmot::ContinuumMechanics::StressMeasures
