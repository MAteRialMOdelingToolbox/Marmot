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
     * @brief Extracts the uniaxial stress-strain tangent from a 3D stiffness matrix.
     *
     * @details Computes the tangent \f$\frac{\partial\sigma}{\partial\varepsilon}\f$
     * in the loading direction (uniaxial) based on the full 3D stiffness matrix \f$\mathbb{C}\f$
     * provided in Voigt notation.
     *
     * @param C 3D stiffness matrix (\f$6\times6\f$) in Voigt notation.
     * @return Uniaxial stress-strain tangent (\f$\frac{\partial\sigma}{\partial\varepsilon}\f$).
     */
    double getUniaxialStressTangent( const Eigen::Ref< const Matrix6d >& C );
  } // namespace UniaxialStress

  namespace PlaneStrain {

    /**
     * @brief Extracts the plane strain stiffness matrix from a 3D stiffness matrix.
     *
     * @details Reduces the full 3D stiffness matrix \f$\mathbb{C}\f$ to the plane strain
     * stiffness matrix corresponding to the in-plane components.
     *
     * @param C 3D stiffness matrix (\f$6\times6\f$) in Voigt notation.
     * @return Plane strain stiffness matrix (\f$3\times3\f$) in Voigt notation.
     */
    Eigen::Matrix3d getPlaneStrainTangent( const Matrix6d& C );

    /**
     * @brief Extracts the plane strain derivative of the stress with respect to the deformation gradient.
     *
     * @details Reduces the derivative of the 3D stress tensor in Voigt notation
     * with respect to the deformation gradient \f$F_{ij}\f$ to the plane strain setting.
     *
     * @param dStressdDeformationGradient3D Derivative of 3D stress with respect to the deformation gradient
     *                                      (\f$6\times3\times3\f$ tensor).
     * @return Reduced derivative of stress with respect to deformation gradient (\f$3\times2\times2\f$ tensor)
     *         in the plane strain setting.
     */
    EigenTensors::Tensor322d reduce3D_dStress_dDeformationGradient(
      const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * @brief Computes the derivative of the 3D strain tensor with respect to the plane strain tensor.
     *
     * @return Matrix (\f$6\times3\f$) representing \f$\frac{\partial \varepsilon_{3D}}{\partial
     * \varepsilon_{plane}}\f$.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStrain();
  } // namespace PlaneStrain

  namespace PlaneStress {

    /**
     * @brief Extracts the plane stress stiffness matrix from a 3D stiffness matrix.
     *
     * @param C 3D stiffness matrix (\f$6\times6\f$) in Voigt notation.
     * @return Plane stress stiffness matrix (\f$3\times3\f$) in Voigt notation.
     */
    Eigen::Matrix3d getPlaneStressTangent( const Matrix6d& C );

    /**
     * @brief Computes the plane stress derivative of the stress with respect to the deformation gradient.
     *
     * @details Reduces the 3D derivative of the stress tensor with respect to the deformation gradient
     * \f$ F_{ij} \f$ to a plane stress setting.
     *
     * @param dStressdDeformationGradient3D Derivative of 3D stress w.r.t. deformation gradient
     *                                      (\f$6 \times 3 \times 3\f$ tensor).
     * @return Reduced derivative of stress w.r.t. deformation gradient (\f$3 \times 2 \times 2\f$ tensor)
     *         for plane stress.
     */
    EigenTensors::Tensor322d compute_dStress_dDeformationGradient(
      const EigenTensors::Tensor633d& dStressdDeformationGradient3D );

    /**
     * @brief Computes the out-of-plane strain correction component for plane stress conditions.
     *
     * @details Determines the correction \f$\varepsilon_{33}^{\text{corr}}\f$ such that the resulting
     * stress satisfies plane stress conditions.
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
     * out-of-plane correction component \f$\varepsilon_{33}^{\text{corr}}\f$ (Voigt notation).
     */
    Marmot::Vector6d planeStressCompensationStrain( const Marmot::Vector6d& elasticStrain, double nu );

    /**
     * @brief Computes the transformation matrix for the plane stress tangent.
     *
     * @details Returns the 6x6 transformation matrix \f$\mathbf{T}\f$ relating arbitrary strain increments
     * to plane stress increments using the material stiffness.
     *
     * @param tangent 6x6 material tangent stiffness matrix (\f$\mathbf{C}\f$).
     * @return 6x6 transformation matrix (\f$\mathbf{T}\f$) for plane stress correction.
     */
    Matrix6d planeStressTangentTransformationMatrix( const Matrix6d& tangent );

    /**
     * @brief Computes the derivative of the 3D strain tensor with respect to the plane stress strain tensor.
     *
     * @param tangent 6x6 material tangent stiffness matrix (\f$\mathbf{C}\f$).
     * @return 6x3 matrix representing \f$\frac{\partial \varepsilon_{3D}}{\partial \varepsilon_{\sigma_{plane}}}\f$.
     */
    Eigen::Matrix< double, 6, 3 > dStrainDStrainPlaneStress( const Matrix6d& tangent );

    /**
     * @brief Computes the derivative of the plane stress tensor with respect to the 3D stress tensor.
     *
     * @return 3x6 matrix representing \f$\frac{\partial \sigma_{plane}}{\partial \sigma_{3D}}\f$.
     */
    Eigen::Matrix< double, 3, 6 > dStressPlaneStressDStress();
  } // namespace PlaneStress
} // namespace Marmot::ContinuumMechanics
