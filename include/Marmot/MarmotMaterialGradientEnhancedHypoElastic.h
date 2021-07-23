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
#include "Marmot/MarmotMaterialGradientEnhancedMechanical.h"

/**
 *
 * Derived abstract base class for gradient-enhanced hypo-elastic materials expressed purely in rate form.
 *
 * In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress tensor
 * \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ the local variable \f$\kappa\f$, the nonlocal variable
 * \f$\bar{\kappa}\f$, internal variables \f$\boldsymbol{x}\f$ and the time \f$ t \f$
 *
 * \f[  \displaystyle \sigRate = f( \sig, \kappa, \epsRate, \dot{\bar{\kappa}}, \boldsymbol{x}, t, ...) \f]
 *
 * In course of numerical time integration, this relation will be formulated incrementally as
 *
 * \f[  \displaystyle \Delta \sig = f ( \sig_n, \kappa_n, \Delta\eps, \Delta\bar{\kappa}, \boldsymbol{x}_n, \Delta t,
 * t_n, ...) \f]
 *
 * with
 *
 * \f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]
 * \f[  \displaystyle \Delta\bar{\kappa} =   \dot{\bar{\kappa}}\, \Delta t \f]
 *
 * and the algorithmic tangents
 *
 * \f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]
 * \f[ \displaystyle \frac{d \sig }{d \bar{\kappa}} =  \frac{d \Delta \sig }{d \Delta \bar{\kappa}} \f]
 * \f[ \displaystyle \frac{d \kappa }{d \eps }=  \frac{d \Delta \kappa }{d \Delta \eps } \f]
 * \f[ \displaystyle \frac{d \kappa }{d \bar{\kappa}}=  \frac{d \Delta \kappa }{d \Delta \bar{\kappa}} =
 * l_{\text{nonlocal}} \f]
 *
 */
class MarmotMaterialGradientEnhancedHypoElastic : public MarmotMaterialGradientEnhancedMechanical {

public:
  using MarmotMaterialGradientEnhancedMechanical::MarmotMaterialGradientEnhancedMechanical;

  /**
   * For a given deformation gradient and a nonlocal variable at the old and the current time,
   * compute the Cauchy stress and the local variable and the algorithmic tangents.
   *
   * @param[in,out]	stress Cauchy stress
   * @param[in,out]	K_local local variable
   * @param[in,out]	nonLocalRadius nonlocal radius representing the tangent \f$\frac{d \Delta \kappa }{d \Delta
   * \bar{\kappa}}\f$
   * @param[in,out]	dStressDDDeformationGradient Derivative of the Cauchy stress tensor with respect to
   * the deformation gradient \f$\boldsymbol{F}\f$
   * @param[in,out]	dK_localDDeformationGradient Derivative of the local variable with respect to
   * the deformation gradient \f$\boldsymbol{F}\f$
   * @param[in,out]	dStressDK Derivative of the Cauchy stress tensor with respect to
   * the nonlocal variable
   * @param[in]	FOld	Deformation gradient at the old (pseudo-)time
   * @param[in]	FNew	Deformation gradient at the current (pseudo-)time
   * @param[in]	KOld	Local variable  at the old (pseudo-)time
   * @param[in]	dK Increment of the local variable
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       stress,
                              double&       KLocal,
                              double&       nonLocalRadius,
                              double*       dStress_dDeformationGradient,
                              double*       dKLocal_dDeformationGradient,
                              double*       dStress_dK,
                              const double* FOld,
                              const double* FNew,
                              const double  KOld,
                              const double  dK,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) override;

  /**
   * For a given deformation gradient and a nonlocal variable at the old and the current time,
   * compute the Cauchy stress and the local variable and the algorithmic tangents.
   *
   * @param[in,out]	stress Cauchy stress
   * @param[in,out]	K_local local variable
   * @param[in,out]	nonLocalRadius nonlocal radius representing the tangent \f$\frac{d \Delta \kappa }{d \Delta
   * \bar{\kappa}}\f$
   * @param[in,out]	dStressDDstrain	Algorithmic tangent representing the derivative of the Cauchy stress tensor with
   * respect to the linearized strain
   * @param[in,out]	dK_localDDStrain Derivative of the local variable with respect to
   * the linearized strain
   * @param[in,out]	dStressDK Derivative of the Cauchy stress tensor with respect to
   * the nonlocal variable
   * @param[in]	KOld	Local variable  at the old (pseudo-)time
   * @param[in]	dK Increment of the local variable
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       stress,
                              double&       K_local,
                              double&       nonLocalRadius,
                              double*       dStressDDStrain,
                              double*       dK_localDDStrain,
                              double*       dStressDK,
                              const double* dStrain,
                              double        KOld,
                              double        dK,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

  /**
   * Plane stress implementation of @ref computeStress.
   *
   * @todo why is dStrain an in-out parameter?
   */
  using MarmotMaterialGradientEnhancedMechanical::computePlaneStress;
  /**
   * Plane stress implementation of @ref computeStress.
   *
   */
  virtual void computePlaneStress( double*       stress2D,
                                   double&       KLocal2D,
                                   double&       nonLocalRadius,
                                   double*       dStress_DStrain2D,
                                   double*       dKLocal_dStrain2D,
                                   double*       dStress_dK2D,
                                   const double* dStrain2D,
                                   double        KOld,
                                   double        dK,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT );
};
