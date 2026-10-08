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
    /**
     * @brief Construct the material and discretize the relaxation function.
     *
     * @param[in] materialProperties properties [E, nu, m, n, nMaxwell, minTau, timeToDays, optional density]
     * @param[in] nMaterialProperties number of properties, at least 7
     * @param[in] materialNumber number of the material
     *
     * @throws std::invalid_argument if fewer than seven properties are given, if there is less than one Maxwell
     * unit, if @p minTau is not positive, or if m < 0 or n <= 0
     */
    LinearViscoElasticWiechert( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /**
     * @brief Compute the stress and the algorithmic tangent for a strain increment.
     *
     * The stress increment is the elastic response with the current effective stiffness, reduced by the stress
     * released by the relaxation of the Maxwell units. The stresses of the Maxwell units are updated afterwards.
     *
     * @param[in,out] state stress and state variables, which hold the stresses of the Maxwell units
     * @param[out] dStressDDStrain algorithmic tangent
     * @param[in] dStrain strain increment in Voigt notation
     * @param[in] timeInfo time and time increment
     */
    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    /**
     * @brief Get the density of the material.
     *
     * @param[in] stateVars state variables, not used
     * @return the density, which is the eighth material property
     *
     * @throws std::runtime_error if the density is not provided
     */
    double getDensity( const double* stateVars ) const override;

  private:
    Wiechert::Properties relaxationTimes;
    Wiechert::Properties elasticModuli;
    double               zerothWiechertStiffness;

    /**
     * @brief Isotropic stiffness tensor for a unit Young's modulus, cached for the instance's lifetime.
     *
     * It depends only on nu, which is fixed at construction, so it is built once here rather than on
     * every computeStress call, as relaxationTimes and elasticModuli already are.
     */
    Marmot::Matrix6d unitStiffness;

    static constexpr int powerLawApproximationOrder = 2;
  };

} // namespace Marmot::Materials
