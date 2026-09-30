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

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::MarmotMeshfreeKernelFunction
   * @brief Abstract interface of a meshfree kernel (window) function @f$ \phi_A(\boldsymbol{x}) @f$.
   *
   * A kernel function is attached to a node (particle) @f$ A @f$ with center @f$ \boldsymbol{x}_A @f$ and has
   * compact support. The reproducing kernel approximation (MarmotMeshfreeReproducingKernelApproximation) builds its
   * shape functions @f$ \Psi_A @f$ by correcting these kernels such that polynomials up to a given order are
   * reproduced exactly; the continuity of @f$ \Psi_A @f$ is inherited from the kernel.
   *
   * All coordinates are passed as raw arrays of length @f$ d @f$ (the spatial dimension of the kernel). The
   * configuration in which the coordinates are given (reference or current) is decided by the caller; the kernel
   * only compares them with its center coordinates.
   */
  class MarmotMeshfreeKernelFunction {

  public:
    /// @brief Default constructor.
    MarmotMeshfreeKernelFunction() = default;

    /// @brief Virtual destructor; kernels are used through base-class pointers.
    virtual ~MarmotMeshfreeKernelFunction() = default;

    /**
     * @brief Evaluate the kernel function @f$ \phi_A(\boldsymbol{x}) @f$.
     * @param[in] coord Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point (length @f$ d @f$).
     * @return The kernel value; zero outside the support.
     */
    virtual double computeKernelFunction( const double* coord ) const = 0;

    /**
     * @brief Evaluate the gradient @f$ \partial \phi_A / \partial x_i @f$ of the kernel function.
     * @param[in]  coord Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point (length @f$ d @f$).
     * @param[out] grad  Gradient of the kernel with respect to @f$ \boldsymbol{x} @f$ (length @f$ d @f$).
     */
    virtual void computeKernelFunctionGradient( const double* coord, double* grad ) const = 0;

    /**
     * @brief Get the axis-aligned bounding box of the support.
     * @param[out] lower Lower corner of the bounding box (length @f$ d @f$).
     * @param[out] upper Upper corner of the bounding box (length @f$ d @f$).
     */
    virtual void getBoundingBox( double* lower, double* upper ) const = 0;

    /**
     * @brief Check whether a point lies within the support of the kernel.
     * @param[in] coord Coordinates of the point (length @f$ d @f$).
     * @return True if the point lies within the support.
     */
    virtual bool isInSupport( const double* coord ) const = 0;

    /**
     * @brief Get the center coordinates @f$ \boldsymbol{x}_A @f$ of the kernel.
     * @return Pointer to the center coordinates (length @f$ d @f$).
     */
    virtual const double* getCenterCoordinates() const = 0;

    /**
     * @brief Move the center of the kernel to a new position.
     * @param[in] coord New center coordinates (length @f$ d @f$).
     */
    virtual void moveTo( const double* coord ) = 0;
  };

} // namespace Marmot::Meshfree
