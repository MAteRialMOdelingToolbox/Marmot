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
#include "Marmot/MarmotMaterialHypoElasticEnergyTangents.h"
#include "Marmot/MarmotStateHelpers.h"

namespace Marmot::Materials {

  /// @brief An implementation of classical J2 plasticity with isotropic hardening.
  class VonMisesModel : public Marmot::Materials::MaterialHypoElastic,
                        public Marmot::Materials::HypoElasticEnergyTangentProvider {

  public:
    /// @brief Numerical constants for the VonMises return-mapping algorithm.
    struct Constants {
      static constexpr double innerNewtonTol        = 1e-12; ///< Convergence tolerance for the inner Newton iteration.
      static constexpr int    nMaxInnerNewtonCycles = 15;    ///< Maximum number of inner Newton iterations.
    };

  public:
    VonMisesModel( const double* materialProperties, const int nMaterialProperties, const int materialLabel );

    void computeStress( state3D&                      state,
                        Marmot::Math::Matrix6d&       dStressDDStrain,
                        const Marmot::Math::Vector6d& dStrain,
                        const timeInfo&               timeInfo ) const override;

    /**
     * @brief computeStress(), which also yields the derivatives of the updated elastic energy density and
     * dissipation wrt the strain increment (Marmot::Materials::HypoElasticEnergyTangentProvider).
     * @param[in,out] state As for computeStress().
     * @param[out] dStress_dStrain As for computeStress().
     * @param[in] dStrain As for computeStress().
     * @param[in] timeInfo As for computeStress().
     * @param[out] dElasticEnergyDensity_dStrain Derivative of the updated elastic energy density.
     * @param[out] dDissipation_dStrain Derivative of the updated dissipation.
     */
    void computeStressWithEnergyTangents( state3D&                      state,
                                          Marmot::Math::Matrix6d&       dStress_dStrain,
                                          const Marmot::Math::Vector6d& dStrain,
                                          const timeInfo&               timeInfo,
                                          Marmot::Math::Vector6d&       dElasticEnergyDensity_dStrain,
                                          Marmot::Math::Vector6d&       dDissipation_dStrain ) const override;

    void computeStressExplicit( state3D&                      state,
                                const Marmot::Math::Vector6d& dStrain,
                                const timeInfo&               timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    /**
     * @brief The implicit stress update of computeStress() and computeStressWithEnergyTangents(), with the check of
     * its results; the energy tangents only if both pointers are given.
     * @param[in,out] state Stress, energies and state variables.
     * @param[out] dStress_dStrain Algorithmic tangent.
     * @param[in] dStrain Strain increment.
     * @param[in] timeInfo Time and time increment.
     * @param[out] dElasticEnergyDensity_dStrain Derivative of the updated elastic energy density, or nullptr.
     * @param[out] dDissipation_dStrain Derivative of the updated dissipation, or nullptr.
     */
    void computeStressImplicit( state3D&                      state,
                                Marmot::Math::Matrix6d&       dStress_dStrain,
                                const Marmot::Math::Vector6d& dStrain,
                                const timeInfo&               timeInfo,
                                Marmot::Math::Vector6d*       dElasticEnergyDensity_dStrain,
                                Marmot::Math::Vector6d*       dDissipation_dStrain ) const;
  };

} // namespace Marmot::Materials
