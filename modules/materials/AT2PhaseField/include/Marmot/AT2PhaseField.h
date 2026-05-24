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
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @brief Simple AT2 phase-field model for fracture simulations.
   *
   * Implements the simple AT2 phase-field fracture model
   * within the general gradient-enhanced hypoelastic framework.
   * See the documentation for details on the theory.
   *
   * **Material properties** (in order):
   *  - \f$ E \f$   : Young's modulus
   *  - \f$ \nu \f$ : Poisson's ratio
   *  - \f$ G_c \f$ : Critical fracture energy
   *  - \f$ l \f$   : Internal length scale
   *  - \f$ \rho \f$ : Mass density (optional, mandatory for dynamic simulations)
   *  - \f$ \eta \f$ : phase-field viscosity (optional, mandatory for dynamic simulations)
   */
  class AT2PhaseField : public MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > {

  public:
    AT2PhaseField( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    void initializeStateLayout()
    {
      stateLayout.add( "maxCrackDrivingForce", 1 );
      stateLayout.add( "strain", 6 );
      stateLayout.finalize();
    }

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    /// @brief Elastic stiffness tensor
    const Marmot::Matrix6d C;
  };

} // namespace Marmot::Materials
