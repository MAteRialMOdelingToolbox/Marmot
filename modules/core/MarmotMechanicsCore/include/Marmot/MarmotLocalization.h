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
#include "Marmot/MarmotTypedefs.h"

/**
 * @file MarmotLocalization.h
 * @brief Localization (strain-softening band) analysis via the acoustic tensor.
 *
 * Functions for checking material instability by evaluating whether the acoustic
 * tensor \f$\mathbf{Q}(\mathbf{n}) = \mathbf{n} \cdot \mathbb{C} \cdot \mathbf{n}\f$
 * loses positive definiteness for some normal vector \f$\mathbf{n}\f$.
 */

namespace Marmot::ContinuumMechanics::LocalizationAnalysis {

  /**
   * @brief Computes the acoustic tensor for a given material tangent and normal vector.
   * @param materialTangent  Fourth-order material tangent in Voigt notation.
   * @param normalVector     Candidate band normal \f$\mathbf{n}\f$.
   * @return Acoustic tensor \f$\mathbf{Q} = \mathbf{n} \cdot \mathbb{C} \cdot \mathbf{n}\f$.
   */
  Marmot::Matrix3d computeAcousticTensor( const Marmot::Matrix6d& materialTangent,
                                          const Marmot::Vector3d& normalVector );

  /**
   * @brief Checks whether the given acoustic tensor indicates localization.
   * @param acousticTensor Acoustic tensor to test.
   * @return @c true if \f$\det(\mathbf{Q}) \le 0\f$ (localization detected).
   */
  bool localizationChecker( const Marmot::Matrix3d& acousticTensor );

  /**
   * @brief Constructs a unit normal vector from spherical angles.
   * @param alpha Polar angle (rad).
   * @param beta  Azimuthal angle (rad).
   * @return Unit normal vector \f$\mathbf{n}\f$.
   */
  Marmot::Vector3d computeNormalVector( double alpha, double beta );

  /**
   * @brief Computes the minimum determinant of the acoustic tensor over all band normals.
   * @param materialTangent Fourth-order material tangent in Voigt notation.
   * @return Minimum determinant of the acoustic tensor.
   */
  double minimumDeterminantAcousticTensor( const Marmot::Matrix6d& materialTangent );

} // namespace Marmot::ContinuumMechanics::LocalizationAnalysis
