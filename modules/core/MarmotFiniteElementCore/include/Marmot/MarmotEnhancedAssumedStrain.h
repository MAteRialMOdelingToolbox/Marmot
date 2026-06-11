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

namespace Marmot {
  /**
   * @brief Enhanced Assumed Strain (EAS) element enrichment utilities.
   *
   * Functions and enumerations for constructing the EAS interpolation matrices
   * used in incompatible-mode and enhanced assumed strain finite element
   * formulations (de Borst, Simo & Rifai variants).
   */
  namespace FiniteElement::EAS {

    /**
     * @brief Supported EAS enrichment types.
     */
    enum EASType {
      DeBorstEAS2,    ///< De Borst 2-parameter EAS
      DeBorstEAS2_P2, ///< De Borst 2-parameter EAS, variant P2
      EAS3,           ///< 3-parameter EAS
      DeBorstEAS6b,   ///< De Borst 6-parameter EAS, variant b
      DeBorstEAS9,    ///< De Borst 9-parameter EAS
      SimoRifaiEAS5,  ///< Simo–Rifai 5-parameter EAS
      SimoRifaiEAS4,  ///< Simo–Rifai 4-parameter EAS
    };

    /**
     * @brief Computes the EAS transformation matrix from the element Jacobian.
     * @param J  Element Jacobian matrix at the reference point.
     * @return   Transformation matrix \f$\mathbf{F}\f$ for mapping EAS modes.
     */
    Eigen::MatrixXd F( const Eigen::MatrixXd& J );

    /**
     * @brief Evaluates the EAS interpolation matrix at a given natural coordinate.
     * @param type EAS enrichment type.
     * @param xi   Natural coordinates of the integration point.
     * @return     EAS interpolation matrix \f$\mathbf{M}\f$.
     */
    Eigen::MatrixXd EASInterpolation( EASType type, const Eigen::VectorXd& xi );

  } // namespace FiniteElement::EAS
} // namespace Marmot
