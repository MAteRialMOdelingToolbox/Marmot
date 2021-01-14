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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot {
    namespace ContinuumMechanics::Kinematics {
        namespace strain {
            Marmot::Vector6d                 GreenLagrange( const Eigen::Matrix3d& F );
            Marmot::EigenTensors::Tensor633d dGreenLagrangedDeformationGradient( const Eigen::Matrix3d& F );
        } // namespace strain

        namespace velocityGradient {
            extern const Eigen::TensorFixedSize< double, Eigen::Sizes< 3, 3, 3, 3 > > dOmega_dVelocityGradient;

            extern const Eigen::TensorFixedSize< double, Eigen::Sizes< 6, 3, 3 > > dStretchingRate_dVelocityGradient;
        } // namespace velocityGradient

        namespace deformationGradient {
            template < int nDim >
            Eigen::Matrix3d make3D( const Eigen::Ref< const Eigen::Matrix< double, nDim, nDim > >& tensor );
        }
    } // namespace ContinuumMechanics::Kinematics
} // namespace Marmot
