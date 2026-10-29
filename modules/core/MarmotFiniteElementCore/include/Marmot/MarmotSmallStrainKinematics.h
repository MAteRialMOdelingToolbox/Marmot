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
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Dense>

namespace Marmot::FiniteElement::SmallStrain {

  /// @brief Strain or stress in (plane) Voigt notation for nDim spatial dimensions (1, 3 or 6 components).
  template < int nDim >
  using Voigt = Eigen::Matrix< double, ContinuumMechanics::VoigtNotation::voigtSizeFromDimension( nDim ), 1 >;

  /**
   * @brief Strain increment from nodal displacement increments, without forming the B-operator.
   * @details The symmetric part of the incremental displacement gradient
   * \f$\partial_\mathbf{x} \Delta\mathbf{u} = \Delta\mathbf{U}\, \partial_\mathbf{x}\mathbf{N}^\mathsf{T}\f$, with
   * the nodal increments \f$\Delta\mathbf{U}\f$ arranged as a nDim x nNodes matrix, in Voigt notation with
   * engineering shear strains. This is \f$\mathbf{B}\, \Delta\mathbf{\qu}\f$ with a sixth (3D) of the data and half
   * of the operations of the dense product.
   * @tparam nDim Number of spatial dimensions (1, 2 or 3).
   * @tparam nNodes Number of element nodes.
   * @param[in] dQU Incremental nodal displacements, node-wise (u1, u2, u3 of node 1, then node 2, ...).
   * @param[in] dNdX Shape function gradients with respect to the physical coordinates (nDim x nNodes).
   * @return The strain increment in (plane) Voigt notation.
   */
  template < int nDim, int nNodes >
  Voigt< nDim > strainIncrement( const Eigen::Ref< const Eigen::Matrix< double, nDim * nNodes, 1 > >& dQU,
                                 const Eigen::Matrix< double, nDim, nNodes >&                         dNdX )
  {
    const Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > > dU( dQU.data() );
    const Eigen::Matrix< double, nDim, nDim >                       dU_dX = dU * dNdX.transpose();

    Voigt< nDim > dE;
    if constexpr ( nDim == 1 )
      dE << dU_dX( 0, 0 );
    else if constexpr ( nDim == 2 )
      dE << dU_dX( 0, 0 ), dU_dX( 1, 1 ), dU_dX( 0, 1 ) + dU_dX( 1, 0 );
    else
      dE << dU_dX( 0, 0 ), dU_dX( 1, 1 ), dU_dX( 2, 2 ), dU_dX( 0, 1 ) + dU_dX( 1, 0 ), dU_dX( 0, 2 ) + dU_dX( 2, 0 ),
        dU_dX( 1, 2 ) + dU_dX( 2, 1 );
    return dE;
  }

  /**
   * @brief Add the internal force of a stress at a quadrature point, \f$\mathbf{B}^\mathsf{T} \boldsymbol{\sigma}\,
   * J_0 w\f$, without forming the B-operator.
   * @details Computed as \f$(\boldsymbol{\sigma} J_0 w)\, \partial_\mathbf{x} \mathbf{N}\f$, a nDim x nNodes matrix
   * of nodal forces (see strainIncrement()).
   * @tparam nDim Number of spatial dimensions (1, 2 or 3).
   * @tparam nNodes Number of element nodes.
   * @param[in,out] fU The internal force vector of the displacement field, node-wise.
   * @param[in] dNdX Shape function gradients with respect to the physical coordinates (nDim x nNodes).
   * @param[in] stress The stress to integrate in (plane) Voigt notation.
   * @param[in] J0xW Integration weight times the Jacobian determinant (times thickness or cross section).
   */
  template < int nDim, int nNodes >
  void addInternalForce( Eigen::Ref< Eigen::Matrix< double, nDim * nNodes, 1 > > fU,
                         const Eigen::Matrix< double, nDim, nNodes >&            dNdX,
                         const Voigt< nDim >&                                    stress,
                         double                                                  J0xW )
  {
    Eigen::Matrix< double, nDim, nDim > sigma;
    if constexpr ( nDim == 1 )
      sigma << stress( 0 );
    else if constexpr ( nDim == 2 )
      sigma << stress( 0 ), stress( 2 ), stress( 2 ), stress( 1 );
    else
      sigma = ContinuumMechanics::VoigtNotation::voigtToStress( stress );

    Eigen::Map< Eigen::Matrix< double, nDim, nNodes > > fUNodal( fU.data() );
    fUNodal.noalias() += ( sigma * J0xW ) * dNdX;
  }

} // namespace Marmot::FiniteElement::SmallStrain
