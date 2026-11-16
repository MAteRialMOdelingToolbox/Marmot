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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

namespace Marmot::Materials {

  /**
   * @brief Multi-surface Rankine (maximum principal stress) plasticity with isotropic hardening/softening.
   *
   * @details Three yield surfaces \f$ f_A = \sigma_A - f_\mathrm{y}(\kappa) \f$, one per principal stress, share
   * the hardening law of VonMisesModel. The return mapping is performed in principal stress space with an
   * active-set strategy; the consistent tangent is assembled from the spectral representation.
   *
   * Material properties: \f$ E, \nu, f_\mathrm{y}^0, H_\mathrm{iso,lin}, \Delta f_\mathrm{y}^{0,\infty}, \delta \f$
   * and optionally the density.
   */
  class RankineSofteningPlasticity : public MaterialHypoElastic {

  public:
    RankineSofteningPlasticity( const double* materialProperties,
                                const int     nMaterialProperties,
                                const int     materialLabel );

    void computeStress( state3D&              state,
                        Math::Matrix6d&       dStressDDStrain,
                        const Math::Vector6d& dStrain,
                        const timeInfo&       timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    /// the stress update of computeStress, without the check of its results
    void updateStress( state3D&              state,
                       Math::Matrix6d&       dStressDDStrain,
                       const Math::Vector6d& dStrain,
                       const timeInfo&       timeInfo ) const;
  };

} // namespace Marmot::Materials
