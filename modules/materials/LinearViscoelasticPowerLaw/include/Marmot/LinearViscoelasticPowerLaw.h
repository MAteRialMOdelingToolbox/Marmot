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

#pragma once
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include <map>
#include <string>

namespace Marmot::Materials {

  /**
   * @brief Implementation of an isotropic linear viscoelastic material
   * using a Power Law compliance function and assuming constant Poisson's ratio
   * generalized for 3D stress states.
   *
   */
  class LinearViscoelasticPowerLaw : public MarmotMaterialHypoElastic {

    /// @brief Young's modulus
    const double& E;
    /**< #E represents the Young's modulus for isotropic linear elasticity.
     * It is a reference variable to #materialProperties[0]. */

    /// @brief Poisson's ratio
    const double& nu;
    /**< #nu represents Poisson's ratio for isotropic linear elasticity.
     * It is a reference variable to #materialProperties[1]. */

    /// @brief power law compliance parameter
    const double& m;
    /**< #m represents the power law compliance parameter.
     * It is a reference variable to #materialProperties[2]. */

    /// @brief power law exponent
    const double& n;
    /**< #n represents the power law exponent.
     * It is a reference variable to #materialProperties[3]. */

    /// @brief number of Kelvin units to approximate the viscoelastic compliance
    const size_t nKelvin;
    /**< #nKelvin represents the number of Kelvin units to approximate the power law function.
     * It is a reference variable to #materialProperties[4]. */

    /// @brief minimal retardation time used in the viscoelastic Kelvin chain
    const double& minTau;
    /**< #minTau represents the minimal retardation time used in the Kelvin chain.
     * It is a reference variable to #materialProperties[5]. */

    /// @brief ratio of simulation time to days
    const double& timeToDays;
    /**< #timeToDays represents the ratio of simulation time to days.
     * It is a reference variable to #materialProperties[6]. */

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearViscoelasticPowerLaw( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( state3D&        state,
                        double*         dStressDDStrain,
                        const double*   dStrain,
                        const timeInfo& timeInfo ) const override;

    void initializeStateLayout() override
    {
      stateLayout.add( "kelvinStateVars", 6 * nKelvin );
      stateLayout.finalize();
    }

  private:
    /// @brief Young's modulus of the #nKelvin Kelvin units
    KelvinChain::Properties elasticModuli;
    /// @brief retardation times of the #nKelvin Kelvin units
    KelvinChain::Properties retardationTimes;
    /// @brief compliance of the zeroth Kelvin unit
    double zerothKelvinChainCompliance;

    static constexpr int powerLawApproximationOrder = 2;
  };
} // namespace Marmot::Materials
