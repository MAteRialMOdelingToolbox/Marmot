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

#include <algorithm>
namespace Marmot::TimeIntegration {

  //! \brief Newmark-Beta time integration
  //! \details
  //! This function implements the Newmark-Beta time integration method for
  //! linear dynamic problems. It updates the velocity and acceleration
  //! vectors based on the displacement increment and the time step.
  //! \param nDim Number of dimensions (2D or 3D)
  //! \param du Displacement increment vector
  //! \param v Velocity vector
  //! \param a Acceleration vector
  //! \param dT Time step
  //! \param newmarkBeta Newmark-Beta parameter (0.25 for linear problems)
  //! \param newmarkGamma Newmark-Gamma parameter (0.5 for linear problems)
  //! \param da_ddu Derivative of acceleration with respect to displacement
  //! \return void
  //! \note The function assumes that the input vectors are of size nDim.
  //! \note The function also assumes that the time step is positive and
  //!       that the Newmark-Beta parameter is not zero.
  //! \note The function uses a small value (1e-16) to avoid division by zero
  template < int nDim >
  void newmarkBetaIntegration( const double* du,
                               double*       v,
                               double*       a,
                               double        dT,
                               double        newmarkBeta,
                               double        newmarkGamma,
                               double*       da_ddu )
  {
    for ( int i = 0; i < nDim; i++ ) {
      dT                    = std::max( dT, 1e-16 );
      const double du_tilde = dT * v[i] + 0.5 * dT * dT * ( ( 1 - 2 * newmarkBeta ) * a[i] );
      const double v_tilde  = v[i] + dT * ( 1 - newmarkGamma ) * a[i];

      a[i]                 = newmarkBeta != 0 ? ( du[i] - du_tilde ) / ( newmarkBeta * dT * dT ) : 0;
      da_ddu[i + i * nDim] = newmarkBeta != 0 ? 1. / ( newmarkBeta * dT * dT ) : 0.0;
      v[i]                 = v_tilde + newmarkGamma * dT * a[i];
    }
  }

} // namespace Marmot::TimeIntegration
