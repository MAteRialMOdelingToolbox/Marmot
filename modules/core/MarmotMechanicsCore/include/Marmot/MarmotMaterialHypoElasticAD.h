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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "autodiff/forward/dual.hpp"

class MarmotMaterialHypoElasticAD : public MarmotMaterialHypoElastic {

public:
  using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

  /**
   * @brief Compute the Cauchy stress tensor \f$\boldsymbol{\sigma}\f$ given an increment of the linearized strain
   * tensor \f$\Delta\boldsymbol{\varepsilon}\f$.
   *
   * Dual numbers are used for both the Cauchy stress tensor and the linearized strain increment, such that the
   * algorithmic tangent operator
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$ can be obtained by means
   * of automatic differentiation.
   *
   * @param[in,out] stress  Cauchy stress tensor
   * @param[in]             dStrain linearized strain increment
   * @param[in]             timeOld Old (pseudo-)time
   * @param[in]             dT (Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]         pNewDT Suggestion for a new time increment
   */
  virtual void computeStressAD( autodiff::dual*       stress,
                                const autodiff::dual* dStrain,
                                const double*         timeOld,
                                const double          dT,
                                double&               pNewDT ) = 0;

  virtual void computeStress( double*       stress,
                              double*       dStressDDStrain,
                              const double* dStrain,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) override;
};
