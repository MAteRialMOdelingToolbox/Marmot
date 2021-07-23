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
 *
 * Derived abstract base class for elastic materials expressed purely in rate form.
 *
 * In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress tensor
 * \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ and the time \f$ t \f$.
 *
 * \f[  \displaystyle \sigRate = f( \sig, \epsRate, t, ...) \f]
 *
 * In course of numerical time integration, this relation will be formulated incrementally as
 *
 * \f[  \displaystyle \Delta \sig = f ( \sig_n, \Delta\eps, \Delta t, t_n, ...) \f]
 *
 * with
 *
 * \f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]
 *
 * and the algorithmic tangent
 *
 * \f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]
 *
 * This formulation is compatible with an Abaqus interface.
 */
class MarmotMaterialHypoElastic : public MarmotMaterialMechanical {

public:
  using MarmotMaterialMechanical::MarmotMaterialMechanical;

  /// Characteristic element length
  double characteristicElementLength;
  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior based on the mesh-adjusted softening
   * modulus.
   *
   * @param[in] length characteristic length; will be assigned to @ref characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  /**
   * For a given deformation gradient at the old and the current time, compute the Cauchy stress and the algorithmic
   * tangent \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{F}^{(n+1)}}\f$.
   *
   * @todo A default implementation is provided.
   *
   * @param[in,out]	stress Cauchy stress
   * @param[in,out]	dSdE	Algorithmic tangent representing the derivative of the Cauchy stress tensor with respect to
   * the deformation gradient \f$\boldsymbol{F}\f$
   * @param[in]	FOld	Deformation gradient at the old (pseudo-)time
   * @param[in]	FNew	Deformation gradient at the current (pseudo-)time
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       stress,
                              double*       dStressDDStrain,
                              const double* FOld,
                              const double* FNew,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) override;

  /**
   * For a given linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$ at the old and the current time,
   * compute the Cauchy stress and the algorithmic tangent
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$.
   *
   * @param[in,out]	stress          Cauchy stress
   * @param[in,out]	dStressDDstrain	Algorithmic tangent representing the derivative of the Cauchy stress tensor with
   * respect to the linearized strain
   * @param[in]	dStrain linearized strain increment
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*       stress,
                              double*       dStressDDStrain,
                              const double* dStrain,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

  /**
   * Plane stress implementation of @ref computeStress.
   */
  using MarmotMaterialMechanical::computePlaneStress;
  virtual void computePlaneStress( double*       stress2D,
                                   double*       dStress_dStrain2D,
                                   const double* dStrain2D,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT );

  /**
   * Uniaxial stress implementation of @ref computeStress.
   */
  using MarmotMaterialMechanical::computeUniaxialStress;
  virtual void computeUniaxialStress( double*       stress1D,
                                      double*       dStress_dStrain1D,
                                      const double* dStrain,
                                      const double* timeOld,
                                      const double  dT,
                                      double&       pNewDT );
};
