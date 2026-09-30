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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot {
  namespace ContinuumMechanics::Kinematics {
    namespace Strain {
      /**
       * @brief Computes the Green-Lagrange strain tensor from the deformation gradient.
       *
       * @param F The deformation gradient tensor as a 3x3 matrix.
       * @return The Green-Lagrange strain tensor represented as a 6D vector in Voigt notation.
       */
      Marmot::Vector6d GreenLagrange( const Eigen::Matrix3d& F );
      /**
       * @brief Computes the derivative of the Green-Lagrange strain tensor with respect to the deformation gradient.
       *
       * @param F The deformation gradient tensor as a 3x3 matrix.
       * @return A \f$6 \times 3 \times 3\f$ tensor representing the derivative of the Green-Lagrange strain tensor.
       */
      Marmot::EigenTensors::Tensor633d dGreenLagrangedDeformationGradient( const Eigen::Matrix3d& F );
    } // namespace Strain

    namespace VelocityGradient {
      /** @brief Computes \f$ \frac{\partial \Omega}{\partial l} \f$ once, for initializing dOmega_dVelocityGradient.
       * @return The 4th-order derivative of the spin tensor with respect to the velocity gradient.
       */
      Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3, 3 > > initializeDOmega_dVelocityGradient();

      /** @brief Computes the derivative of the stretching rate with respect to the velocity gradient once, for
       * initializing dStretchingRate_dVelocityGradient.
       * @return The \f$ 6 \times 3 \times 3 \f$ derivative, stretching rate in engineering Voigt notation.
       */
      Eigen::TensorFixedSize< double, Eigen::Sizes< 6, 3, 3 > > initializeDStretchingRate_dVelocityGradient();

      /**
       * @brief A 4th-order tensor representing the derivative of the spin tensor \f$ \Omega \f$
       *        with respect to the velocity gradient tensor.
       */
      inline const Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3, 3 > >
        dOmega_dVelocityGradient = initializeDOmega_dVelocityGradient();

      /** @brief Tensor representing the derivative of the stretching rate with respect to the velocity gradient.*/
      inline const Eigen::TensorFixedSize< double, Eigen::Sizes< 6, 3, 3 > >
        dStretchingRate_dVelocityGradient = initializeDStretchingRate_dVelocityGradient();
    } // namespace VelocityGradient

    namespace DeformationGradient {
      /**
       * @brief Embeds a lower-dimensional square tensor into a 3x3 tensor.
       *
       * This function takes an \f$n \times n\f$ tensor (with \f$n=1,2,3\f$)
       * and returns a 3x3 tensor:
       * - For \f$n=1\f$: the scalar is placed in the (0,0) entry of a 3x3 identity matrix.
       * - For \f$n=2\f$: the 2x2 tensor is placed in the top-left block of a 3x3 identity matrix.
       * - For \f$n=3\f$: the tensor is returned unchanged.
       *
       * @tparam nDim Dimension of the input tensor (1, 2, or 3).
       * @param tensor An \f$n \times n\f$ tensor represented as Eigen::Matrix<double,nDim,nDim>.
       * @return A 3x3 tensor embedding the input according to the rules above.
       */
      template < int nDim >
      Eigen::Matrix3d make3D( const Eigen::Ref< const Eigen::Matrix< double, nDim, nDim > >& tensor );
    } // namespace DeformationGradient
  }   // namespace ContinuumMechanics::Kinematics
} // namespace Marmot
