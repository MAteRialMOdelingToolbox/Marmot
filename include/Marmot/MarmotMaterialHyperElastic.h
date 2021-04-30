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
#include "Marmot/MarmotMaterialMechanical.h"

/**
 * Derived abstract base class for _simple_, purely hyperelastic materials to be used for finite elements based on the
 * total lagrangian kinematic description (TL elements). The second Piola - Kirchhoff stress tensor \f$ S \f$ will be
 * derived by
 *
 * \f[ \displaystyle S = \frac{\partial f(\boldsymbol{E},t )}{\partial \boldsymbol{E}} \f]
 *
 * with the Green - Lagrange strain tensor \f$ \boldsymbol{E} \f$
 *
 * \f[
 *   \displaystyle E  = \frac{1}{2}\,\left(\boldsymbol{F}^T\cdot \boldsymbol{F} - \boldsymbol{I} \right)
 * \f]
 *
 * as work conjugated measure and the variable \f$ \boldsymbol{F} \f$ denoting the deformation gradient.
 * The algorithmic tangent will be calculated by
 *
 * \f[
 *   \displaystyle \frac{d \boldsymbol{S}}{d \boldsymbol{E}}
 * \f]
 */

class MarmotMaterialHyperElastic : public MarmotMaterialMechanical {

public:
  using MarmotMaterialMechanical::MarmotMaterialMechanical;
  /**
   *
   * For a given deformation gradient at the old and the current time, compute the 2nd Piola-Kirchhoff stress and the
   * algorithmic tangent \f$\frac{\partial\boldsymbol{S}^{(n+1)}}{\partial\boldsymbol{E}^{(n+1)}}\f$.
   *
   * @todo A default implementation is provided.
   *
   * @param[in,out]	S	2nd Piola-Kirchhoff stress
   * @param[in,out]	dSdE	Algorithmic tangent representing the derivative of the 2nd Piola-Kirchhoff stress tensor with
   * respect to the Green-Lagrange strain tensor \f$\boldsymbol{E}\f$.
   * @param[in]	FOld	Deformation gradient at the old (pseudo-)time
   * @param[in]	FNew	Deformation gradient at the current (pseudo-)time
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       S,    // PK2
                              double*       dSdE, // d PK2 d GL_E
                              const double* FOld,
                              const double* FNew,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) override;

  /**
   * For a given Green-Lagrange strain, compute the 2nd Piola-Kirchhoff stress and the algorithmic tangent
   * \f$\frac{\partial\boldsymbol{S}^{(n+1)}}{\partial\boldsymbol{E}^{(n+1)}}\f$.
   *
   * @todo Should we use function overloading in this case and simple also use computeStress for the function name?
   *
   * @param[in,out]	S	2nd Piola-Kirchhoff stress
   * @param[in,out]	dSdE	Algorithmic tangent representing the derivative of the 2nd Piola-Kirchhoff stress tensor with
   * respect to the Green-Lagrange strain tensor \f$\boldsymbol{E}\f$.
   * @param[in]	deltaE Green-Lagrange strain increment
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStressPK2( double*       S,
                                 double*       dSdE,
                                 const double* E,
                                 const double* timeOld,
                                 const double  dT,
                                 double&       pNewDT ) = 0;

  /**
   * Plane stress implementation of @ref computeStressPK2.
   */
  virtual void computePlaneStressPK2( double*       S,
                                      double*       dSdE,
                                      double*       E,
                                      const double* timeOld,
                                      const double  dT,
                                      double&       pNewDT );

  /**
   * Uniaxial stress implementation of @ref computeStressPK2.
   */
  virtual void computeUniaxialStressPK2( double*       S,
                                         double*       dSdE,
                                         double*       E,
                                         const double* timeOld,
                                         const double  dT,
                                         double&       pNewDT );
};
