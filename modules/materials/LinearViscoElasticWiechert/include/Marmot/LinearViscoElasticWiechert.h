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
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
   * A power-law relaxation function \f$\Psi(t)=m t^{-n}\f$ is approximated
   * by logarithmically spaced Maxwell branches. `E` is the equilibrium
   * Young's modulus.
   *
   * Material properties are ordered as
   * `[E, nu, m, n, nMaxwell, minTau, timeToDays, optional density]`.
   */
  class LinearViscoElasticWiechert : public MarmotMaterialHypoElastic {

    /// \brief Young's modulus
    const double& E;

    /// \brief Poisson's ratio
    const double& nu;

    /// \brief scaling factor of the power-law relaxation function
    const double& m;

    /// \brief exponent of the power-law relaxation function
    const double& n;

    /// \brief number of Maxwell units used to approximate the relaxation function
    const size_t nMaxwell;

    /// \brief minimal relaxation time used in the viscoelastic Maxwell chain
    const double& minTau;

    /// \brief ratio of simulation time to days
    const double& timeToDays;

  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    LinearViscoElasticWiechert( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void initializeStateLayout()
    {
      stateLayout.add( "maxwellStateVars", 6 * nMaxwell );
      stateLayout.finalize();
    }

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    /// @brief relaxation times of the #nMaxwell Maxwell units
    Wiechert::Properties relaxationTimes;
    /// @brief Young's modulus of the #nMaxwell Maxwell units
    Wiechert::Properties elasticModuli;
    /// @brief stiffness of the zeroth Wiechert unit
    double zerothWiechertStiffness;

    static constexpr int powerLawApproximationOrder = 2;
  };

} // namespace Marmot::Materials
