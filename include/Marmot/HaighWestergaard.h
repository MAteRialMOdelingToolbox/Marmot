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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotTypedefs.h"

namespace Marmot {
  namespace ContinuumMechanics::HaighWestergaard {

    /**
     * Aggregate of the Haigh-Westergaard coordinates (invariants) \f$\xi\f$,
     * \f$\rho\f$, \f$\theta\f$.
     */
    struct HaighWestergaardCoordinates {

      /// Hydrostatic component \f$\xi\f$
      double xi;
      ///  Deviatoric radius \f$\rho\f$
      double rho;
      /// Lode angle \f$\theta\f$ specified in radian
      double theta;
    };

    /**
     * Computes the stress coordinates in the Haigh-Westergaard space.
     *
     * \note The stress coordinates are computed from the invariants \f$I_1,\,J_2,\,J_3\f$ of the stress tensor as
     * follows:
     *
     * \f[\xi=I_1/\sqrt{3}\f]
     * \f[\rho=\sqrt{2\,J_2}\f]
     * \f[\theta=\frac{1}{3}\,\arccos{\left(\frac{3\,\sqrt{3}}{2}\,\frac{J_3}{\sqrt{J_2^3}}\right)}\f]
     *
     * @param stress Stress tensor \f$\sig\f$ given in \ref voigtnotation "Voigt notation".
     */
    HaighWestergaardCoordinates haighWestergaard( const Marmot::Vector6d& stress );
    /**
     * Computes the strain coordinates in the Haigh-Westergaard space.
     *
     * \note The computation is equal to @ref haighWestergaard by replacing the stress invariants with the strain
     * invariants.
     *
     * @param strain Strain tensor \f$\eps\f$ given in \ref voignotation "Voigt notation".
     */
    HaighWestergaardCoordinates haighWestergaardFromStrain( const Marmot::Vector6d& strain );

  } // namespace ContinuumMechanics::HaighWestergaard
} // namespace Marmot
