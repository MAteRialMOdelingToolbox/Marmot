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
#include "Marmot/MarmotPortability.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include <string>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::BondSlipLaw
   * @brief Abstract base class for bond-slip laws, i.e., the constitutive relation between the slip of a
   * reinforcement bar relative to its host continuum and the bond stress acting on the bar surface.
   *
   * A bond-slip law operates in the local frame of the bar: the first component is the tangential slip along the bar
   * axis, the second and third components are the relative displacements normal to the bar axis. In two dimensions
   * the third component is zero. The returned bond stress is a traction per unit bar surface; the element integrates it
   * over the bar perimeter and length.
   *
   * Bond-slip laws are pluggable: they live as modules in the materials category and register themselves in the
   * Marmot::Registry::BondSlipLawFactory, from which bond elements create them by name.
   */
  class MARMOT_API BondSlipLaw {

  protected:
    const double* materialProperties;  ///< Pointer to the array of material properties
    const int     nMaterialProperties; ///< Number of material properties

  public:
    const int materialNumber; ///< Integer identifier for this material instance

    /**
     * @brief Construct the bond-slip law with a given set of properties and an identifier.
     * @param[in] materialProperties_  Pointer to the array of material properties.
     * @param[in] nMaterialProperties_ Number of entries in @p materialProperties_.
     * @param[in] materialNumber_      Integer identifying this material instance.
     */
    BondSlipLaw( const double* materialProperties_, int nMaterialProperties_, int materialNumber_ )
      : materialProperties( materialProperties_ ),
        nMaterialProperties( nMaterialProperties_ ),
        materialNumber( materialNumber_ )
    {
    }

    /// Default destructor
    virtual ~BondSlipLaw() = default;

    /// Layout of the state variables
    Marmot::Utilities::StateLayoutDynamic stateLayout;

    /// Bond state at an integration point of the bar: values at the beginning of the increment on entry, updated
    /// values on exit
    struct State {
      Marmot::Math::Vector3d bondStress;           ///< [tangential bond stress, normal stress 1, normal stress 2]
      double                 elasticEnergyDensity; ///< elastic energy per unit bar surface
      double                 dissipation;          ///< dissipated energy per unit bar surface
      double*                stateVars;            ///< pointer to the state variables
    };

    /// Time at the beginning of the increment and size of the time increment
    struct TimeIncrement {
      double time; ///< time at the beginning of the increment
      double dT;   ///< size of the time increment
    };

    /**
     * @brief Compute the bond stress and its algorithmic tangent.
     * @param[inout] state  values at the beginning of the increment on entry, updated values on exit
     * @param[out] dBondStress_dSlip algorithmic tangent \f$\frac{\partial \boldsymbol{\tau}}{\partial
     * \boldsymbol{s}}\f$
     * @param[in] slip   total slip at the end of the increment, local bar frame [tangential, normal 1, normal 2]
     * @param[in] dSlip  slip increment, local bar frame
     * @param[in] timeIncrement time information
     */
    virtual void computeBondStress( State&                        state,
                                    Marmot::Math::Matrix3d&       dBondStress_dSlip,
                                    const Marmot::Math::Vector3d& slip,
                                    const Marmot::Math::Vector3d& dSlip,
                                    const TimeIncrement&          timeIncrement ) const = 0;

    /**
     * @brief Get a view to the state variables.
     * @param stateName Name of the state variable
     * @param stateVars Pointer to the state variable array
     * @return Marmot::Utilities::StateView to access the state variable
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, double* stateVars ) const
    {
      return stateLayout.getStateView( stateVars, stateName );
    }

    /**
     * @brief Get the total number of required state variables.
     * @return Total number of required state variables
     */
    int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

    /**
     * @brief Initialize the state variables.
     * @param stateVars Pointer to the state variable array
     * @param nStateVars Number of state variables
     *
     * @note The default implementation initializes all state variables to zero.
     */
    virtual void initializeYourself( double* stateVars, int nStateVars )
    {
      for ( int i = 0; i < nStateVars; ++i )
        stateVars[i] = 0.0;
    }
  };

} // namespace Marmot::Materials
