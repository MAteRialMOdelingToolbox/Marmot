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

    /**
     * @brief Computes the Kirchhoff stress from a principal Kirchhoff stress given in the
     * eigenbasis of the right Cauchy-Green tensor.
     *
     * @details For an isotropic hyperelastic material formulated in principal stretches
     * \f$\lambda_i\f$, the second Piola-Kirchhoff stress shares the eigenbasis \f$\boldsymbol Q\f$
     * of \f$\boldsymbol C=\boldsymbol F^T\boldsymbol F\f$ with the principal Kirchhoff stress via
     * \f$S_i = \tau_i/\lambda_i^2\f$. This helper assembles
     * \f$\boldsymbol S = \boldsymbol Q\,\mathrm{diag}(S_i)\,\boldsymbol Q^T\f$ and pushes it
     * forward to the Kirchhoff stress \f$\boldsymbol\tau=\boldsymbol F\boldsymbol S\boldsymbol
     * F^T\f$, factoring out the pattern shared by every principal-stretch-based finite-strain
     * material (e.g. Marmot::Materials::Ogden and its incompressible Neo-Hooke/Mooney-Rivlin
     * specializations).
     *
     * @tparam T Scalar type, e.g. double, autodiff::dual.
     * @param principalKirchhoffStress Principal Kirchhoff stresses \f$\tau_i\f$.
     * @param principalStretches Principal stretches \f$\lambda_i\f$ (same ordering/eigenbasis).
     * @param Q Eigenvector matrix of \f$\boldsymbol C\f$ (e.g. from
     * DeformationMeasures::principalStretches()).
     * @param F Deformation gradient.
     * @return Kirchhoff stress \f$\boldsymbol\tau\f$.
     */
    template < typename T >
    Tensor33t< T > KirchhoffStressFromPrincipalKirchhoffStress( const Tensor3t< T >&  principalKirchhoffStress,
                                                                const Tensor3t< T >&  principalStretches,
                                                                const Tensor33t< T >& Q,
                                                                const Tensor33t< T >& F )
    {
      Tensor33t< T > PK2_principal( 0. );
      for ( int i = 0; i < 3; ++i )
        PK2_principal( i, i ) = principalKirchhoffStress( i ) / ( principalStretches( i ) * principalStretches( i ) );

      const Tensor33t< T > PK2 = Q % PK2_principal % transpose( Q );

      return KirchhoffStressFromPK2( PK2, F );
    }

    /**
     * @brief Converts a Biot stress given in the eigenbasis of the right stretch tensor to the
     * second Piola-Kirchhoff stress.
     *
     * @details For an isotropic hyperelastic (or hyperelastic-viscoelastic) material formulated in
     * terms of the Biot stress \f$\boldsymbol S_0\f$ work-conjugate to the right stretch tensor
     * \f$\boldsymbol U=\boldsymbol Q\,\mathrm{diag}(\lambda_i)\,\boldsymbol Q^T\f$, the second
     * Piola-Kirchhoff stress follows from rotating \f$\boldsymbol S_0\f$ into the eigenbasis
     * \f$\boldsymbol Q\f$ of \f$\boldsymbol U\f$ (and of \f$\boldsymbol C=\boldsymbol
     * F^T\boldsymbol F\f$), converting via \f$S_{ij} = 2\,(S_0)^{\rm rot}_{ij}/(\lambda_i+
     * \lambda_j)\f$, and rotating back:
     * \f[
     *   \boldsymbol S = \boldsymbol Q\,\mathrm{diag\_pair}\!\left(\frac{2\,(\boldsymbol
     *   Q^T\boldsymbol S_0\boldsymbol Q)_{ij}}{\lambda_i+\lambda_j}\right)\boldsymbol Q^T .
     * \f]
     * Unlike KirchhoffStressFromPrincipalKirchhoffStress(), \f$\boldsymbol S_0\f$ need not be
     * diagonal in the \f$\boldsymbol Q\f$ basis (e.g. after a history-dependent/viscoelastic
     * update), hence the pairwise \f$\lambda_i+\lambda_j\f$ conversion rather than a per-component
     * \f$\lambda_i^2\f$ one - the two helpers are not interchangeable.
     *
     * @tparam T Scalar type, e.g. double, autodiff::dual.
     * @param biotStress Biot stress \f$\boldsymbol S_0\f$, in the reference (unrotated) basis.
     * @param principalStretches Principal stretches \f$\lambda_i\f$ of \f$\boldsymbol U\f$ (same
     * ordering/eigenbasis as \p Q).
     * @param Q Eigenvector matrix of \f$\boldsymbol U\f$ (e.g. from
     * DeformationMeasures::principalStretches()).
     * @return Second Piola-Kirchhoff stress \f$\boldsymbol S\f$.
     */
    template < typename T >
    Tensor33t< T > PK2FromRotatedBiotStress( const Tensor33t< T >& biotStress,
                                             const Tensor3t< T >&  principalStretches,
                                             const Tensor33t< T >& Q )
    {
      const Tensor33t< T > biotStressRotated = transpose( Q ) % biotStress % Q;

      Tensor33t< T > PK2_rotated( 0. );
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          PK2_rotated( i, j ) = 2.0 * biotStressRotated( i, j ) / ( principalStretches( i ) + principalStretches( j ) );

      return Q % PK2_rotated % transpose( Q );
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
