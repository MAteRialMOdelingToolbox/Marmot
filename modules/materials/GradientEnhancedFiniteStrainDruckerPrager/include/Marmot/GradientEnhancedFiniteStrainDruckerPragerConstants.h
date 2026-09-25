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
 * Thomas Mader thomas.mader@boku.ac.at
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

namespace Marmot::Materials {

  namespace GradientEnhancedFiniteStrainDruckerPragerConstants {
    /// absolute tolerance of the return-map residual (log strains are dimensionless, the yield
    /// function is scaled by the cohesive strength)
    const double innerNewtonTol        = 1e-12;
    const int    nMaxInnerNewtonCycles = 25;
    /// forward-difference step of the algorithmic tangents, relative to the unit entries of F
    const double tangentPerturbation = 1e-7;
    /// forward-difference step of the tangents w.r.t. the nonlocal field
    const double tangentPerturbationNonLocal = 1e-8;
    /// the deviatoric Mandel stress counts as vanished (apex) below this fraction of the cohesion
    const double apexTol = 1e-12;
  } // namespace GradientEnhancedFiniteStrainDruckerPragerConstants

} // namespace Marmot::Materials
