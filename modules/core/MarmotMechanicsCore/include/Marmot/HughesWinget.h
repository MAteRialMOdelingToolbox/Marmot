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
#include "Marmot/MarmotVoigt.h"

/**
 * @file HughesWinget.h
 * @brief Hughes–Winget objective stress-rate integrator.
 *
 * Provides the @c HughesWinget class that, given old and new deformation
 * gradients, computes the objective strain increment, the incremental rotation
 * tensor, and the tangent operator contributions required by hypoelastic
 * constitutive updates.
 */

namespace Marmot::NumericalAlgorithms {
  /**
   * @brief Hughes–Winget objective integration algorithm.
   *
   * Given the deformation gradients at the beginning and end of an increment,
   * this class computes:
   * - the symmetric strain increment in Voigt notation,
   * - the incremental rotation matrix,
   * - a rotated stress tensor, and
   * - consistent tangent contributions for mixed deformation-gradient/stress
   *   quantities.
   */
  class HughesWinget {
  public:
    /**
     * @brief Available integration formulations.
     */
    enum Formulation {
      AbaqusLike ///< Formulation consistent with the Abaqus UMAT conventions.
    };

    /**
     * @brief Constructs a HughesWinget integrator from two successive deformation gradients.
     * @param FOld       Deformation gradient at the start of the increment.
     * @param FNew       Deformation gradient at the end of the increment.
     * @param formulation Integration formulation to use.
     */
    HughesWinget( const Eigen::Matrix3d& FOld, const Eigen::Matrix3d& FNew, Formulation formulation );

    /// @brief Returns the objective symmetric strain increment in Voigt notation.
    Marmot::Vector6d getStrainIncrement();
    /// @brief Returns the incremental rotation matrix \f$\Delta\mathbf{R}\f$.
    Eigen::Matrix3d getRotationIncrement();
    /**
     * @brief Rotates a Voigt stress (or strain) tensor with the incremental rotation.
     * @param tensor Input tensor in Voigt notation.
     * @return Rotated tensor in Voigt notation.
     */
    Marmot::Vector6d rotateTensor( const Marmot::Vector6d& tensor );

    /**
     * @brief Compute the consistent tangent contribution \f$\frac{d\boldsymbol{S}}{d\mathbf{F}}\f$.
     *
     * Given the current Cauchy stress, the inverse deformation gradient, and the
     * Cauchy-stress/strain tangent, returns the third-order tensor mapping variations
     * of the deformation gradient to variations of the Cauchy stress.
     *
     * @param stress        Current Cauchy stress in Voigt notation.
     * @param FInv          Inverse of the deformation gradient at the end of the increment.
     * @param dChauchyDEps  Algorithmic Cauchy-stress/strain tangent (6×6).
     * @return              The 6×3×3 tangent tensor \f$\frac{d\boldsymbol{S}}{d\mathbf{F}}\f$.
     */
    Marmot::EigenTensors::Tensor633d compute_dS_dF( const Marmot::Vector6d& stress,
                                                    const Eigen::Matrix3d&  FInv,
                                                    const Marmot::Matrix6d& dChauchyDEps );

    /**
     * @brief Compute the consistent tangent \f$\frac{d\phi}{d\mathbf{F}}\f$ for a scalar quantity.
     *
     * Given the inverse deformation gradient and the derivative of a scalar (e.g., an
     * internal variable) with respect to the strain in Voigt notation, returns the
     * corresponding 3×3 matrix mapping variations of the deformation gradient.
     *
     * @param FInv         Inverse of the deformation gradient at the end of the increment.
     * @param dScalarDEps  Derivative of the scalar quantity with respect to the strain, in Voigt notation.
     * @return             The 3×3 sensitivity matrix \f$\frac{d\phi}{d\mathbf{F}}\f$.
     */
    Eigen::Matrix3d compute_dScalar_dF( const Eigen::Matrix3d& FInv, const Marmot::Vector6d& dScalarDEps );

  private:
    Formulation      theFormulation;
    Eigen::Matrix3d  l;
    Eigen::Matrix3d  dOmega;
    Eigen::Matrix3d  dR;
    Marmot::Vector6d dEps;
  };
} // namespace Marmot::NumericalAlgorithms
