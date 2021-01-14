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

namespace Marmot::NumericalAlgorithms {
    class HughesWinget {
      public:
        enum Formulation {
            AbaqusLike

        };

        HughesWinget( const Eigen::Matrix3d& FOld, const Eigen::Matrix3d& FNew, Formulation formulation );

        Marmot::Vector6d getStrainIncrement();
        Eigen::Matrix3d  getRotationIncrement();
        Marmot::Vector6d rotateTensor( const Marmot::Vector6d& tensor );

        Marmot::EigenTensors::Tensor633d compute_dS_dF( const Marmot::Vector6d& stress,
                                                        const Eigen::Matrix3d&  FInv,
                                                        const Marmot::Matrix6d& dChauchyDEps );
        Eigen::Matrix3d compute_dScalar_dF( const Eigen::Matrix3d& FInv, const Marmot::Vector6d& dScalarDEps );

      private:
        Formulation      theFormulation;
        Eigen::Matrix3d  l;
        Eigen::Matrix3d  dOmega;
        Eigen::Matrix3d  dR;
        Marmot::Vector6d dEps;
    };
} // namespace Marmot::NumericalAlgorithms
