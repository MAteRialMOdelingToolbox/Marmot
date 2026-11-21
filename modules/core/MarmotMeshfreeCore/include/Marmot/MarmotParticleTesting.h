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
#pragma once

#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Core>
#include <string>
#include <vector>

/**
 * @brief Helpers for the tests of particles.
 */
namespace Marmot::Testing {

  /**
   * @brief Check the bound of @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$ a particle gives
   * (Particle::computeBoundOfAbsoluteStiffnessTimes()) against its stiffness @f$ \boldsymbol{K} @f$: the bound never
   * lies below @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$, and the largest entry of each column (what the
   * Collatz-Wielandt bound takes) is at most @p maxOverestimate times the exact one.
   * @param name The name of the particle, for the messages.
   * @param particle The particle, at the state of @p K.
   * @param stateVars Its state variables, reset after the bound as a host does.
   * @param K The stiffness of the particle at a zero increment, column-major, nDof x nDof.
   * @param W The nonnegative weights, nDof x nColumns.
   * @param maxOverestimate The largest admissible ratio of the largest entries of a column.
   * @throws std::runtime_error (throwExceptionOnFailure()) if the bound fails either check.
   */
  inline void checkAbsoluteStiffnessBoundAgainstStiffness( const std::string&                           name,
                                                           Marmot::Discretizations::Meshfree::Particle& particle,
                                                           std::vector< double >&                       stateVars,
                                                           const Eigen::MatrixXd&                       K,
                                                           const Eigen::MatrixXd&                       W,
                                                           double maxOverestimate )
  {
    const int       nColumns = static_cast< int >( W.cols() );
    const auto      backup   = stateVars;
    Eigen::MatrixXd bound    = Eigen::MatrixXd::Zero( W.rows(), nColumns );
    particle.computeBoundOfAbsoluteStiffnessTimes( W.data(), nColumns, bound.data(), 1.0, 1.0 );
    stateVars = backup;

    const Eigen::MatrixXd exact = K.cwiseAbs() * W;
    for ( int c = 0; c < nColumns; c++ ) {
      const double largest = exact.col( c ).maxCoeff();
      throwExceptionOnFailure( ( bound.col( c ) - exact.col( c ) ).minCoeff() >= -1e-12 * largest,
                               Marmot::Utilities::MakeString()
                                 << name << ": the bound of |K| W lies below it in column " << c );
      const double overestimate = bound.col( c ).maxCoeff() / largest;
      throwExceptionOnFailure( overestimate <= maxOverestimate,
                               Marmot::Utilities::MakeString()
                                 << name << ": the bound of |K| W overestimates its largest entry of column " << c
                                 << " by " << overestimate );
    }
  }

} // namespace Marmot::Testing
