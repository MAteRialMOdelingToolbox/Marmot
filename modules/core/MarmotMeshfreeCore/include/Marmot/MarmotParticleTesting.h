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

#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleStiffness.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Core>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Helpers for the tests of particles.
 */
namespace Marmot::Testing {

  /**
   * @brief Check the bound of @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$ of a particle
   * (Marmot::Discretizations::Meshfree::addBoundOfAbsoluteStiffnessTimes()) against its stiffness @f$ \boldsymbol{K}
   * @f$: the bound never lies below @f$ |\boldsymbol{K}|\, \boldsymbol{W} @f$, and the largest entry of each column
   * (what the Collatz-Wielandt bound takes) is at most @p maxOverestimate times the exact one.
   * @param name The name of the particle, for the messages.
   * @param particle The particle, at the state of @p K.
   * @param K The stiffness of the particle at a zero increment, column-major, nDof x nDof.
   * @param W The nonnegative weights, nDof x nColumns.
   * @param maxOverestimate The largest admissible ratio of the largest entries of a column.
   * @throws std::runtime_error (throwExceptionOnFailure()) if the bound fails either check.
   */
  inline void checkAbsoluteStiffnessBoundAgainstStiffness( const std::string&                           name,
                                                           Marmot::Discretizations::Meshfree::Particle& particle,
                                                           const Eigen::MatrixXd&                       K,
                                                           const Eigen::MatrixXd&                       W,
                                                           double maxOverestimate )
  {
    const int       nColumns = static_cast< int >( W.cols() );
    Eigen::MatrixXd bound    = Eigen::MatrixXd::Zero( W.rows(), nColumns );
    Marmot::Discretizations::Meshfree::addBoundOfAbsoluteStiffnessTimes( particle,
                                                                         W.data(),
                                                                         nColumns,
                                                                         bound.data(),
                                                                         1.0,
                                                                         1.0 );

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

  /**
   * @brief Check that a particle that describes its stiffness (Particle::describeStiffness()) is left as it was found,
   * in the middle of an increment: advanced by @p dQ and not yet accepted, so that its state holds the increment and
   * its material point keeps the response of it, both of which the zero increment of the description would overwrite.
   * The state variables, the bound of two descriptions in a row (addBoundOfAbsoluteStiffnessTimes()) and, for a
   * particle with the explicit path, the internal force before and after a description, which reads the response the
   * material point keeps, must be equal to the last bit.
   * @param name The name of the particle, for the messages.
   * @param particle The particle, with kernel functions assigned.
   * @param stateVars Its state variables.
   * @param dQ The increment of its nodes, node-wise.
   * @param explicitPath Whether the particle has the explicit path: it is advanced by updatePhysicsExplicit() and its
   * explicit internal force is compared; else it is advanced by computePhysicsKernels().
   * @throws std::runtime_error (throwExceptionOnFailure()) if the description changes the particle.
   */
  inline void checkDescribingTheStiffnessLeavesTheParticleUnchanged(
    const std::string&                           name,
    Marmot::Discretizations::Meshfree::Particle& particle,
    const std::vector< double >&                 stateVars,
    const Eigen::VectorXd&                       dQ,
    bool                                         explicitPath )
  {
    const auto nDof = dQ.size();
    if ( explicitPath )
      particle.updatePhysicsExplicit( dQ.data(), 1.0, 1.0 );
    else {
      Eigen::VectorXd force   = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd tangent = Eigen::MatrixXd::Zero( nDof, nDof );
      particle.computePhysicsKernels( dQ.data(), force.data(), tangent.data(), 1.0, 1.0 );
    }

    const auto forceAndBound = [&]() {
      Eigen::VectorXd force = Eigen::VectorXd::Zero( nDof );
      if ( explicitPath )
        particle.computePhysicsKernelsExplicit( force.data() );
      const Eigen::VectorXd W     = Eigen::VectorXd::Ones( nDof );
      Eigen::VectorXd       bound = Eigen::VectorXd::Zero( nDof );
      Marmot::Discretizations::Meshfree::addBoundOfAbsoluteStiffnessTimes( particle,
                                                                           W.data(),
                                                                           1,
                                                                           bound.data(),
                                                                           1.0,
                                                                           1.0 );
      return std::make_pair( force, bound );
    };

    const std::vector< double > stateBefore   = stateVars;
    const auto [forceBefore, boundBefore]     = forceAndBound();
    const std::vector< double > stateAfterOne = stateVars;
    const auto [forceAfterOne, boundAfterOne] = forceAndBound();
    const std::vector< double > stateAfterTwo = stateVars;

    throwExceptionOnFailure( stateAfterOne == stateBefore && stateAfterTwo == stateBefore,
                             name + ": describing the stiffness changes the state variables" );
    throwExceptionOnFailure( forceAfterOne == forceBefore,
                             name + ": describing the stiffness changes the internal force that follows it" );
    throwExceptionOnFailure( boundAfterOne == boundBefore,
                             name + ": a second description of the stiffness differs from the first" );
  }

} // namespace Marmot::Testing
