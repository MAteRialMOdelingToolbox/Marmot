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

#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {
  /**
   * @brief Shrinkage strain computation helpers for the B4 model.
   */
  namespace Shrinkage {
    /**
     * @brief B4-specific shrinkage strain increment computation.
     */
    namespace B4 {

      /**
       * @brief Computes the total shrinkage strain increment over one time step.
       *
       * Evaluates the combined autogenous and drying shrinkage strain increment
       * according to the B4 model.
       *
       * @param tStartDays                     Age at the start of the increment (days).
       * @param dTDays                         Length of the time increment (days).
       * @param ultimateAutogenousShrinkageStrain Ultimate autogenous shrinkage strain
       * \f$\eps^\mathrm{shr,au}_\infty\f$.
       * @param autogenousShrinkageHalfTime    Autogenous shrinkage half-time \f$\tau_\mathrm{shr,au}\f$ (days).
       * @param alpha                          Exponent parameter \f$\alpha\f$ for autogenous shrinkage.
       * @param rt                             Drying shrinkage half-time \f$\tau_\mathrm{shr,d}\f$ (days).
       * @param ultimateDryingShrinkageStrain  Ultimate drying shrinkage strain at zero humidity
       * \f$\eps^\mathrm{shr,d}_\infty\f$.
       * @param dryingShrinkageHalfTime        Alternate drying shrinkage half-time parameter (days).
       * @param kHum                           Humidity-dependent factor \f$k_\mathrm{h}\f$.
       * @param dryingStart                    Age at the start of drying \f$t_0\f$ (days).
       * @return Voigt-notation shrinkage strain increment vector.
       */
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
    } // namespace B4
  }   // namespace Shrinkage
} // namespace Marmot::Materials
