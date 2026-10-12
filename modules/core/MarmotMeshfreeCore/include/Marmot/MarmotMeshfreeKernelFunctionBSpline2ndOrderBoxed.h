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

#include "Marmot/MarmotMeshfreeKernelFunction.h"
#include <cmath>

namespace Marmot::Meshfree {

  /**
   * @class Marmot::Meshfree::MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed
   * @brief Boxed tensor-product kernel from 2nd order (quadratic) B-splines.
   *
   * The kernel of node @f$ A @f$ with center @f$ \boldsymbol{x}_A @f$ and support radius @f$ a @f$ is the
   * tensor product of one-dimensional 2nd order (quadratic) B-splines @f$ w @f$ in each coordinate direction,
   * @f[
   *   \phi_A(\boldsymbol{x}) = \prod_{i=1}^{d} w\left( x_i - x_{A,i} \right),
   * @f]
   * with
   * @f[
   *   w(r) = \begin{cases} 1 - 2 z^2, & 0 \le z \le \tfrac{1}{2}, \\
   *                        2 (1 - z)^2, & \tfrac{1}{2} < z \le 1, \\
   *                        0, & z > 1, \end{cases}
   *   \qquad z = \frac{|r|}{a}.
   * @f]
   * @f$ w @f$ is the quadratic B-spline on the knots @f$ \{-a, -a/2, a/2, a\} @f$, scaled to @f$ w(0) = 1 @f$; it
   * is @f$ C^1 @f$-continuous, so the kernel and the resulting reproducing kernel shape functions are
   * @f$ C^1 @f$.
   *
   * The support is the open box @f$ |x_i - x_{A,i}| < a @f$, @f$ i = 1, \dots, d @f$ (hence "boxed"); its
   * bounding box is @f$ [\boldsymbol{x}_A - a, \boldsymbol{x}_A + a] @f$. The support radius @f$ a @f$ is
   * the half edge length of the box, identical in all directions, and the kernel is not normalized (the reproducing
   * kernel correction removes any constant scaling of the kernel).
   *
   * @note The kernel does not own its center coordinates: it stores the pointer passed to the constructor, and
   * moveTo() writes the new center through this pointer. The pointed-to storage must outlive the kernel.
   */
  class MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed : public MarmotMeshfreeKernelFunction {

  private:
    double*      _centerCoord;   ///< non-owning pointer to the center coordinates @f$ \boldsymbol{x}_A @f$
    const double _supportRadius; ///< support radius @f$ a @f$ (half edge length of the support box)
    const int    _dim;           ///< spatial dimension @f$ d @f$

  public:
    /**
     * @brief Construct the kernel.
     * @param[in] centerCoord   Pointer to the center coordinates @f$ \boldsymbol{x}_A @f$ (length @f$ d @f$); the
     *                          kernel keeps (and moveTo() modifies) this storage, it is not copied.
     * @param[in] dim           Spatial dimension @f$ d @f$.
     * @param[in] supportRadius Support radius @f$ a @f$.
     */
    MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed( double* centerCoord, int dim, double supportRadius );

    /// @brief Default destructor.
    virtual ~MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed() = default;

    /**
     * @brief Evaluate the tensor-product kernel @f$ \phi_A(\boldsymbol{x}) = \prod_i w(x_i - x_{A,i}) @f$.
     * @param[in] coord Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point (length @f$ d @f$).
     * @return The kernel value; zero outside the support box.
     */
    double computeKernelFunction( const double* coord ) const override;

    /**
     * @brief Evaluate the kernel gradient by the product rule.
     *
     * @f$ \partial \phi_A / \partial x_i = w'(x_i - x_{A,i}) \prod_{j \neq i} w(x_j - x_{A,j}) @f$.
     *
     * @param[in]  coord Coordinates @f$ \boldsymbol{x} @f$ of the evaluation point (length @f$ d @f$).
     * @param[out] grad  Gradient @f$ \partial \phi_A / \partial x_i @f$ (length @f$ d @f$).
     */
    void computeKernelFunctionGradient( const double* coord, double* grad ) const override;

    /**
     * @brief Evaluate the one-dimensional B-spline @f$ w(r) @f$.
     * @param[in] coord_minus_center Signed distance @f$ r = x_i - x_{A,i} @f$ in one coordinate direction.
     * @return @f$ w(r) @f$.
     */
    double computeBSpline2ndOrder( double coord_minus_center ) const;

    /**
     * @brief Evaluate the derivative @f$ \mathrm{d}w / \mathrm{d}r @f$ of the one-dimensional B-spline.
     *
     * @f[
     *   \frac{\mathrm{d}w}{\mathrm{d}r} = \frac{\operatorname{sign}(r)}{a} \begin{cases} -4 z, & z \le \tfrac{1}{2},
     *   \\ -4 + 4 z, & \tfrac{1}{2} < z \le 1, \\ 0, & z > 1. \end{cases}
     * @f]
     * @param[in] coord_minus_center Signed distance @f$ r = x_i - x_{A,i} @f$ in one coordinate direction.
     * @return @f$ \mathrm{d}w / \mathrm{d}r @f$.
     */
    double computeBSpline2ndOrderGradient( double coord_minus_center ) const;

    /**
     * @brief Get the center coordinates.
     * @return The (non-owning) pointer passed to the constructor.
     */
    const double* getCenterCoordinates() const override;

    /**
     * @brief Move the kernel center by overwriting the referenced center coordinates.
     * @param[in] coord New center coordinates (length @f$ d @f$).
     */
    void moveTo( const double* coord ) override;

    /**
     * @brief Check whether a point lies within the support, i.e., whether @f$ \phi_A(\boldsymbol{x}) > 0 @f$.
     * @param[in] coord Coordinates of the point (length @f$ d @f$).
     * @return True if the kernel is positive at the point.
     */
    bool isInSupport( const double* coord ) const override;

    /**
     * @brief Get the bounding box @f$ [\boldsymbol{x}_A - a, \boldsymbol{x}_A + a] @f$ of the support.
     * @param[out] min Lower corner @f$ x_{A,i} - a @f$ (length @f$ d @f$).
     * @param[out] max Upper corner @f$ x_{A,i} + a @f$ (length @f$ d @f$).
     */
    void getBoundingBox( double* min, double* max ) const override;
  };

} // namespace Marmot::Meshfree
