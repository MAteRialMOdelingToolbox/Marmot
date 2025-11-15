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
#include "Marmot/MarmotMaterial.h"
#include "Marmot/MarmotTypedefs.h"

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
class MarmotMaterialHypoElastic : public MarmotMaterial {

public:
  using MarmotMaterial::MarmotMaterial;

  /// Structure to hold the material state at a material point in 3D
  struct state3D {
    Marmot::Vector6d stress;       ///> Cauchy stress tensor in Voigt notation
    double           strainEnergy; ///> Strain energy density
    double*          stateVars;    ///> Pointer to array of state variables
  };

  // Structure to hold the material state at a material point for 2D plane stress
  struct state2D {
    Marmot::Vector3d stress;       ///> 2D Cauchy stress tensor in Voigt notation
    double           strainEnergy; ///> Strain energy density
    double*          stateVars;    ///> Pointer to array of state variables
  };

  // Structure to hold the material state at a material point for 1D uniaxial stress
  struct state1D {
    double  stress;       ///> 1D Cauchy stress
    double  strainEnergy; ///> Strain energy density
    double* stateVars;    ///> Pointer to array of state variables
  };

  struct timeInfo {
    double time; ///> Current (pseudo-)time
    double dT;   ///> (Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
  };

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
   * For a given linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$ at the old and the current time,
   * compute the Cauchy stress and the algorithmic tangent
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$.
   *
   * @param[in,out]	state  A state3D instance carrying stress, strain energy, and state variables
   * @param[in,out]	dStressDDstrain	Algorithmic tangent representing the derivative of the Cauchy stress tensor with
   * respect to the linearized strain
   * @param[in]	dStrain linearized strain increment
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   */
  virtual void computeStress( state3D&        state,
                              double*         dStressDDStrain,
                              const double*   dStrain,
                              const timeInfo& timeInfo ) = 0;

  /**
   * Plane stress implementation of @ref computeStress.
   */
  virtual void computePlaneStress( state2D&        stress2D,
                                   double*         dStress_dStrain2D,
                                   const double*   dStrain2D,
                                   const timeInfo& timeInfo );

  /**
   * Uniaxial stress implementation of @ref computeStress.
   */
  virtual void computeUniaxialStress( state1D&        stress1D,
                                      double*         dStress_dStrain1D,
                                      const double*   dStrain,
                                      const timeInfo& timeInfo );
};
