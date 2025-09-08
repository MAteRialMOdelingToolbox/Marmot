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

namespace Marmot::ContinuumMechanics {

  namespace UniaxialStress {
    /**
     * Extract the uniaxial stress-strain tangent \f$\frac{\partial\sigma}{\partial\varepsilon}\f$
     * from a given three-dimensional stiffness matrix.
     *
     * @param C 3D stiffness matrix \f$\mathbb{C}\f$ given in \ref voigtnotation "Voigt notation".
     */
    double getUniaxialStressTangent( const Eigen::Ref< const Matrix6d >& C );
  } // namespace UniaxialStress

  namespace PlaneStrain {

    /**
     * Extract the plane strain stiffness matrix from a given three-dimensional stiffness matrix.
     *
     * @param C 3D stiffness matrix \f$\mathbb{C}\f$ given in \ref voigtnotation "Voigt notation".
     */
    Eigen::Matrix3d getPlaneStrainTangent( const Matrix6d& C );

    /**
     * Extract the plane strain derivitive of the stress in Voigt notation with respect to the
     * deformation gradient \f$F_{ij}\f$ from the corresponding derivative in a 3d setting.
     */
    EigenTensors::Tensor322d reduce3D_dStress_dDeformationGradient(
      const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * Compute the derivative of the three-dimensional strain tensor with respect to the plane
     * strain tensor.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStrain();
  } // namespace PlaneStrain

  namespace PlaneStress {

    /**
     * @brief Extract the plane stress stiffness matrix from a given three-dimensional stiffness matrix.
     *
     * @param C 3D stiffness matrix \f$ \mathbb{C} \f$ given in \ref voigtnotation "Voigt notation".
     */
    Eigen::Matrix3d getPlaneStressTangent( const Matrix6d& C );

    /**
     * @brief Extract the plane stress derivitive of the stress in Voigt notation with respect to the deformation
     * gradient \f$ F_{ij} \f$ from the corresponding derivative in a 3d setting.
     */
    EigenTensors::Tensor322d compute_dStress_dDeformationGradient(
      const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * @brief Computes the out-of-plane strain correction component for plane stress conditions.
     *
     * This function calculates the compensation strain \f$ \varepsilon_{33}^\text{corr} \f$ for a given elastic strain
     * \f$ \boldsymbol{\varepsilon}^{\text{el}} \f$, so that the resulting stress satisfies plane stress:
     * \f[
     *     \sigma^\text{plane} = \mathbb{C} : (\boldsymbol{\varepsilon}^{\text{el}} +
     * \boldsymbol{\varepsilon}^{\text{corr}})
     * \f]
     *
     * @param elasticStrain The in-plane elastic strain vector \f$ \boldsymbol{\varepsilon}^{\text{el}} \f$ (6
     * components, Voigt notation).
     * @param nu Poisson's ratio of the material.
     * @return A 6-component strain vector \f$ \boldsymbol{\varepsilon}^{\text{corr}} \f$ containing the computed
     * out-of-plane correction component (Voigt notation).
     */
    Marmot::Vector6d planeStressCompensationStrain( const Marmot::Vector6d& elasticStrain, double nu );

    /**
     * @brief Computes the transformation matrix for the plane stress tangent.
     *
     * This function returns the 6x6 transformation matrix \f$\mathbf{T}\f$ that relates
     * an arbitrary strain increment to a plane stress increment via the material stiffness:
     * \f[
     *     \Delta \boldsymbol{\sigma}^{\text{plane}} = \mathbb{C} : (\mathbf{T} : \Delta
     *     \boldsymbol{\varepsilon}^{\text{arbitrary}})
     * \f]
     *
     * @param tangent The 6x6 material tangent stiffness matrix \f$\mathbf{C}\f$.
     * @return The 6x6 transformation matrix \f$\mathbf{T}\f$ for plane stress correction.
     */
    Matrix6d planeStressTangentTransformationMatrix( const Matrix6d& tangent );

    /**
     * @brief Compute the derivative of the three-dimensional strain tensor with respect to the strain
     * tensor in a plane stress setting.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStress( const Matrix6d& tangent );

    /**
     * Compute the derivative of the plane stress tensor with respect to the three-dimensional stress tensor.
     */
    Eigen::Matrix< double, 3, 6 > dStressPlaneStressDStress();
  } // namespace PlaneStress
} // namespace Marmot::ContinuumMechanics
