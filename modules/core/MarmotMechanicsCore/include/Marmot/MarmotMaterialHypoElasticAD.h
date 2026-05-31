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

/**
 * @brief Derived class of MarmotMaterialHypoElastic providing automatic differentiation support
 *        for computing the algorithmic tangent via dual numbers.
 */
class MarmotMaterialHypoElasticAD : public MarmotMaterialHypoElastic {

public:
  /// @brief Inherits the base-class constructor; see MarmotMaterialHypoElastic::MarmotMaterialHypoElastic().
  using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

  /// @brief Structure holding the material state for 3D using dual numbers (for automatic differentiation).
  struct state3DAD {
    Marmot::Vector6dual stress;               ///< Cauchy stress tensor in Voigt notation
    double              elasticEnergyDensity; ///< Elastic strain energy density
    double              dissipation;          ///< Dissipation
    double*             stateVars;            ///< Pointer to array of state variables
    /**
     * @brief Default constructor for initializing state3DAD
     * Initializes stress to zero, energy and dissipation to zero, and stateVars to nullptr.
     */
    state3DAD()
      : stress( Marmot::Vector6dual::Zero() ), elasticEnergyDensity( 0.0 ), dissipation( 0.0 ), stateVars( nullptr ){};
    /**
     * @brief Constructor for initializing state3DAD
     *@param state A state3D object containing the initial values for stress, energy density, dissipation, and state
     *variables.
     */
    state3DAD( const state3D& state )
      : stress( state.stress ),
        elasticEnergyDensity( state.elasticEnergyDensity ),
        dissipation( state.dissipation ),
        stateVars( state.stateVars ){};
  };
  /**
   * @brief Compute the Cauchy stress tensor \f$\boldsymbol{\sigma}\f$ given an increment of the linearized strain
   * tensor \f$\Delta\boldsymbol{\varepsilon}\f$.
   *
   * Dual numbers are used for both the Cauchy stress tensor and the linearized strain increment, such that the
   * algorithmic tangent operator
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$ can be obtained by means
   * of automatic differentiation.
   *
   * @param[in,out] state  State carrying the dual Cauchy stress, strain energy, and state variables
   * @param[in]             dStrain linearized strain increment
   * @param[in]             timeInfo Old (pseudo-)time
   */
  virtual void computeStressAD( state3DAD&                 state,
                                const Marmot::Vector6dual& dStrain,
                                const timeInfo&            timeInfo ) const = 0;

  virtual void computeStress( state3D&                state,
                              Marmot::Matrix6d&       dStressDDStrain,
                              const Marmot::Vector6d& dStrain,
                              const timeInfo&         timeInfo ) const override;
};
