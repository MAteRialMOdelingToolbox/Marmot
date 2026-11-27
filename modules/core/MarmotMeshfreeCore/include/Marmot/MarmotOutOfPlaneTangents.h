/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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

namespace Marmot::Discretizations::MaterialPoints {

  /**
   * @brief The out-of-plane parts of the tangent of a 2D material point with respect to the increment of the
   * deformation gradient, which the hoop terms of an axisymmetric particle need (the out-of-plane stretch
   * @f$ F_{33} @f$ is the hoop stretch there).
   *
   * @details With @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$, @f$ \partial F_{iI} / \partial
   * \Delta F_{jJ} = \delta_{ij}\, F_{n,JI} @f$, and a block-diagonal 2D state (no out-of-plane shear),
   * @f[
   *   \frac{\partial \tau_{ij}}{\partial \Delta F_{33}} = \frac{\partial \tau_{ij}}{\partial F_{3I}}\, F_{n,3I},
   *   \qquad
   *   \frac{\partial \tau_{33}}{\partial \Delta F_{kL}} = \frac{\partial \tau_{33}}{\partial F_{kI}}\, F_{n,LI} .
   * @f]
   */
  struct OutOfPlaneTangents {
    /// @f$ \partial \tau_{ij} / \partial \Delta F_{33} @f$, in-plane @f$ ij @f$
    Fastor::Tensor< double, 2, 2 > dTau_dDeltaF33 = Fastor::Tensor< double, 2, 2 >( 0.0 );
    /// @f$ \partial \tau_{33} / \partial \Delta F_{kL} @f$, in-plane @f$ kL @f$
    Fastor::Tensor< double, 2, 2 > dTau33_dDeltaF = Fastor::Tensor< double, 2, 2 >( 0.0 );
    /// @f$ \partial \tau_{33} / \partial \Delta F_{33} @f$
    double dTau33_dDeltaF33 = 0.0;
  };

  /**
   * @brief The out-of-plane parts of the tangent of a 2D material point, see OutOfPlaneTangents.
   * @tparam DeformationGradient A 3 x 3 tensor or a map of one (e.g. on the state variables).
   * @param[in] dTau_dF The 3D tangent @f$ \partial \tau_{ij} / \partial F_{kL} @f$ of the material.
   * @param[in] F_n The 3D deformation gradient of the last accepted state.
   * @return The out-of-plane parts.
   */
  template < typename DeformationGradient >
  OutOfPlaneTangents computeOutOfPlaneTangents( const Marmot::Math::Tensors::FastorTensors::Tensor3333d& dTau_dF,
                                                const DeformationGradient&                               F_n )
  {
    OutOfPlaneTangents t;
    for ( int I = 0; I < 3; I++ ) {
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < 2; j++ ) {
          t.dTau_dDeltaF33( i, j ) += dTau_dF( i, j, 2, I ) * F_n( 2, I );
          t.dTau33_dDeltaF( i, j ) += dTau_dF( 2, 2, i, I ) * F_n( j, I );
        }
      t.dTau33_dDeltaF33 += dTau_dF( 2, 2, 2, I ) * F_n( 2, I );
    }
    return t;
  }

} // namespace Marmot::Discretizations::MaterialPoints
