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

#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {
  namespace Shrinkage {
    namespace B4 {

      Marmot::Vector6d computeShrinkageStrainIncrement( const double tStartDays,
                                                        const double dTDays,
                                                        const double ultimateAutogenousShrinkageStrain,
                                                        const double autogenousShrinkageHalfTime,
                                                        const double alpha,
                                                        const double rt,
                                                        const double ultimateDryingShrinkageStrain,
                                                        const double dryingShrinkageHalfTime,
                                                        const double kHum,
                                                        const double dryingStart );
    }
  } // namespace Shrinkage
} // namespace Marmot::Materials
