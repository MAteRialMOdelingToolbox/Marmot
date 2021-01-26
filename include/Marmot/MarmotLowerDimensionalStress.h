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

namespace Marmot::ContinuumMechanics {

    namespace UniaxialStress {
        double getUniaxialStressTangent( const Eigen::Ref< const Matrix6d >& C );
    }

    namespace PlaneStrain {

        Eigen::Matrix3d getPlaneStrainTangent( const Matrix6d& C );

        EigenTensors::Tensor322d dStressdDeformationGradient(
            const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

        Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStrain();
    } // namespace PlaneStrain

    namespace PlaneStress {

        Eigen::Matrix3d getPlaneStressTangent( const Matrix6d& C );

        EigenTensors::Tensor322d dStressdDeformationGradient(
            const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

        /*compute E33 for a given elastic strain, to compute the compensation for
         * planeStress = Cel : (elasticStrain + compensationStrain) */
        Marmot::Vector6d planeStressCompensationStrain( const Marmot::Vector6d& elasticStrain, double nu );
        /* Returns the transformation Matrix T which fullfills
         * planeStressIncrement = C : (T : arbitraryStrainIncrement) */
        Matrix6d                      planeStressTangentTransformationMatrix( const Matrix6d& tangent );
        Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStress( const Matrix6d& tangent );
        Eigen::Matrix< double, 3, 6 > dStressPlaneStressDStress();
    } // namespace PlaneStress

} // namespace Marmot::ContinuumMechanics
