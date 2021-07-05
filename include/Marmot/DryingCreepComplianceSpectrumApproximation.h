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
#include "Marmot/Solidification.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  namespace SpectrumApproximation {

    namespace DryingCreepB3 {

      enum class ApproximationOrder : int { firstOrder = 1 };

      void computeApproximation( Eigen::Ref< Solidification::KelvinProperties > kelvinElasticModuli,
                                 Eigen::Ref< Solidification::KelvinProperties > kelvinRetardationTimes,
                                 const double                                   hEnv,
                                 const double                                   xi,
                                 const double                                   xiZero,
                                 const enum ApproximationOrder                  approximationOrder );

      Solidification::KelvinProperties generateRetardationTimes( double minTime, int numUnits );

      dualDouble phi( dualDouble xi, double b, double xiZero );
      
      dualDouble T( dualDouble eta );

      dualDouble S( dualDouble xi );

      dualDouble f( dualDouble xi, double b );
      
    } // namespace DryingCreepB3
  }   // namespace SpectrumApproximation
} // namespace Marmot::Materials
