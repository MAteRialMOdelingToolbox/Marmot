/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ ___ ___   ___ | |_
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

#pragma once

#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotWiechert.h"

#include <cstddef>

namespace Marmot::Materials {

  /**
   * @brief Isotropic linear viscoelastic material using a generalized Maxwell chain.
   *
   * A power-law relaxation function Psi(t)=m t^(-n) is approximated by logarithmically spaced Maxwell branches.
   * `E` is the equilibrium Young's modulus.
   *
   * Material properties are ordered as [E, nu, m, n, nMaxwell, minTau, timeToDays, optional density].
   */
  class LinearViscoElasticWiechert : public MarmotMaterialHypoElastic {

    const double& E;
    const double& nu;
    const double& m;
    const double& n;
    const size_t  nMaxwell;
    const double& minTau;
    const double& timeToDays;

  public:
    LinearViscoElasticWiechert( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    Wiechert::Properties relaxationTimes;
    Wiechert::Properties elasticModuli;
    double               zerothWiechertStiffness;

    static constexpr int powerLawApproximationOrder = 2;
  };

} // namespace Marmot::Materials
