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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMaterialHypoElasticAD.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/dual.hpp"
#include <Eigen/src/Core/Map.h>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot;

namespace Marmot::Materials {
  /**
   * @brief Implementation of a isotropic J2-plasticity  material
   * for 3D stress states using automatic differentiation.
   */
  class ADVonMises : public MarmotMaterialHypoElasticAD {
  public:
    using MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD;

    /// @brief Young’s modulus.
    const double& E;
    /// @brief Poisson’s ratio.
    const double& nu;
    /// @brief Initial yield stress.
    const double& yieldStress;
    /// @brief First hardening parameter.
    const double& HLin;
    /// @brief Second hardening parameter.
    const double& deltaYieldStress;
    /// @brief Third hardening parameter.
    const double& delta;
    /// @brief Shear modulus.
    const double G;

    ADVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    void computeStressAD( state3DAD& state, const autodiff::dual* dStrain, const timeInfo& timeInfo ) override;

    class ADVonMisesModelStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kappa", .length = 1 },
      } );

      /// @brief Hardening variable.
      double& kappa;

      ADVonMisesModelStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), kappa( find( "kappa" ) ){};
    };
    std::unique_ptr< ADVonMisesModelStateVarManager > managedStateVars;

    int getNumberOfRequiredStateVars() override { return ADVonMisesModelStateVarManager::layout.nRequiredStateVars; }

    void assignStateVars( double* stateVars, int nStateVars ) override;

    StateView getStateView( const std::string& result ) override;

    /**
     * @brief Hardening function.
     * @tparam[in] kappa Hardening variable.
     * @returns Current yield stress.
     */
    template < typename T >
    T fy( T kappa_ )
    {
      const T res = yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - exp( -delta * kappa_ ) );
      return res;
    }

    /**
     * @brief Yield function.
     * @tparam[in] rho Deviatoric radius (\f$\rho = ||s||\f$).
     * @param[in] kappa Hardening variable.
     * @returns Value of the yield function.
     */
    template < typename T >
    T f( const T rho_, const double kappa_ )
    {
      return rho_ - Constants::sqrt2_3 * fy( kappa_ );
    }

    /**
     * @brief Objective function for stress update algorithm.
     * @tparam[in] rhoTrial Deviatoric radius.
     * @param[in] kappa Current value of the hardening variable.
     * @tparam[in] deltaKappa Increment of the hardening variable.
     * @returns Value of the yield function.
     */
    template < typename T >
    T g( const T rhoTrial, const double kappa, const T deltaKappa )
    {
      const T kappa_ = kappa + deltaKappa;
      return rhoTrial - Constants::sqrt6 * G * deltaKappa - Constants::sqrt2_3 * fy( kappa_ );
    }
  };
} // namespace Marmot::Materials
