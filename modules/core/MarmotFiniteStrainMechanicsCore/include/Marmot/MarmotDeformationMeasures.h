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

    /** @brief Computes the right Cauchy-Green tensor \f$\boldsymbol{C}\f$.
     * @tparam T Scalar type (e.g., float, double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return The right Cauchy-Green tensor
     *
     * The right Cauchy-Green deformation tensor is defined as:
     * \f[
     *   \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F}
     * \f]
     * or in index notation:
     * \f[
     *   C_{IJ} = F_{iI} F_{iJ}
     * \f]
     */
    template < typename T >
    Tensor33t< T > rightCauchyGreen( const Tensor33t< T > F )
    {
      Tensor33t< T > C = einsum< iI, iJ >( F, F );
      return C;
    }

    /** @brief Computes the left Cauchy-Green tensor \f$\boldsymbol{b}\f$.
     * @tparam T Scalar type (e.g., float, double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return The left Cauchy-Green tensor b
     *
     * The left Cauchy-Green deformation tensor is defined as:
     * \f[
     *   \boldsymbol{b} = \boldsymbol{F} \boldsymbol{F}^T
     * \f]
     * or in index notation:
     * \f[
     *   b_{ij} = F_{iI} F_{jI}
     * \f]
     */
    template < typename T >
    Tensor33t< T > leftCauchyGreen( const Tensor33t< T > F )
    {
      Tensor33t< T > b = einsum< iI, JI >( F, F );
      return b;
    }

    namespace FirstOrderDerived {

      /** @brief Computes the right Cauchy-Green tensor \f$\boldsymbol{C}\f$ and its derivative with respect to
       * \f$\boldsymbol{F}\f$.
       * @tparam T Scalar type (e.g., float, double, autodiff::dual)
       * @param F Deformation gradient tensor
       * @returns A pair containing the right Cauchy-Green tensor C and its derivative with respect to F
       *
       * The right Cauchy-Green deformation tensor is defined as:
       * \f[
       *   \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F}
       * \f]
       * or in index notation:
       * \f[
       *   C_{IJ} = F_{iI} F_{iJ}
       * \f]
       * The derivative of \f$\boldsymbol{C}\f$ with respect to \f$\boldsymbol{F}\f$ is given by:
       * \f[
       *  \frac{\partial C_{IJ}}{\partial F_{kK}} = \delta_{IK} F_{kJ} + F_{kI} \delta_{JK}
       * \f]
       */
      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > rightCauchyGreen( const Tensor33t< T > F )
      {
        Tensor33t< T > C = einsum< iI, iJ >( F, F );

        const Tensor3333t< T > dC_dF = einsum< IK, kJ, to_IJkK >( Spatial3D::I, F ) +
                                       einsum< kI, JK, to_IJkK >( F, Spatial3D::I );

        return { C, dC_dF };
      }

      /** @brief Computes the left Cauchy-Green tensor \f$\boldsymbol{b}\f$ and its derivative with respect to
       * \f$\boldsymbol{F}\f$.
       * @tparam T Scalar type (e.g., float, double, autodiff::dual)
       * @param F Deformation gradient tensor
       * @return A pair containing the left Cauchy-Green tensor b and its derivative with respect to F
       *
       * The left Cauchy-Green deformation tensor is defined as:
       * \f[
       *   \boldsymbol{b} = \boldsymbol{F} \boldsymbol{F}^T
       * \f]
       * or in index notation:
       * \f[
       *   b_{ij} = F_{iI} F_{jI}
       * \f]
       * The derivative of \f$\boldsymbol{b}\f$ with respect to \f$\boldsymbol{F}\f$ is given by:
       * \f[
       *  \frac{\partial b_{ij}}{\partial F_{kK}} = \delta_{ik} F_{jK} + F_{iK} \delta_{jk}
       * \f]
       */
      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > leftCauchyGreen( const Tensor33t< T > F )
      {
        Tensor33t< T > b = einsum< iI, JI >( F, F );

        const Tensor3333t< T > db_dF = einsum< ik, jK, to_ijkK >( Spatial3D::I, F ) +
                                       einsum< iK, jk, to_ijkK >( F, Spatial3D::I );

        return { b, db_dF };
      }

    } // namespace FirstOrderDerived

  }   // namespace DeformationMeasures
} // namespace Marmot::ContinuumMechanics
