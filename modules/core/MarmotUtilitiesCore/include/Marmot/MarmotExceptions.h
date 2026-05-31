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
#include <stdexcept>
#include <string>

namespace Marmot {

  /** @class StressUpdateFailed
   *  @brief Exception thrown when a stress update (e.g., return-mapping iteration) fails to converge
   *         inside a material model.
   */
  class StressUpdateFailed : public std::runtime_error {
  public:
    /**
     * @brief Constructor for StressUpdateFailed exception.
     * @param message A descriptive error message explaining the reason for the exception.
     */
    explicit StressUpdateFailed( const std::string& message ) : std::runtime_error( message ) {}
  };

  /** @class SolverConvergenceFailed
   *  @brief Exception thrown when the Newton-Raphson iteration in a material-point solver fails to
   *         converge (maximum number of iterations reached or NaN encountered).
   */
  class SolverConvergenceFailed : public std::runtime_error {
  public:
    /**
     * @brief Constructor for SolverConvergenceFailed exception.
     * @param message A descriptive error message explaining the reason for the exception.
     */
    explicit SolverConvergenceFailed( const std::string& message ) : std::runtime_error( message ) {}
  };

  /** @class SolverTimestepExhausted
   *  @brief Exception thrown when the material-point solver cannot reduce the time step below the
   *         prescribed minimum.
   */
  class SolverTimestepExhausted : public std::runtime_error {
  public:
    /**
     * @brief Constructor for SolverTimestepExhausted exception.
     * @param message A descriptive error message explaining the reason for the exception.
     */
    explicit SolverTimestepExhausted( const std::string& message ) : std::runtime_error( message ) {}
  };

  /** @class SolverIncrementsExhausted
   *  @brief Exception thrown when the material-point solver exhausts the maximum number of
   *         increments without completing the step.
   */
  class SolverIncrementsExhausted : public std::runtime_error {
  public:
    /**
     * @brief Constructor for SolverIncrementsExhausted exception.
     * @param message A descriptive error message explaining the reason for the exception.
     */
    explicit SolverIncrementsExhausted( const std::string& message ) : std::runtime_error( message ) {}
  };

} // namespace Marmot
