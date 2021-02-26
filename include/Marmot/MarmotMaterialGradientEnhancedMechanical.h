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
#include "Marmot/MarmotMaterial.h"

class MarmotMaterialGradientEnhancedMechanical : public MarmotMaterial {

  /*
     Abstract basic class for Mechanical materials.
     'Mechanical' is meant in the 'most general sense', i.e., any material which describes a mechanical (cauchy)
     stress - deformation relationship (e.g, hyperelastic, hypoelastic, elasto-plastic, visco-elastic materials)

     σ = f (σ, dxdX, t, .. ),

     formulated incrementally as σ_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, .. )

     Algorithmic tangent: dσdF = d σ_np d (dxdX_np)

      Format:

      / d σ_11 d F_00,   d σ_11 d F_10,   d σ_11 d F_20,   d σ_11 d F_01, \
      | d σ_22 d F_00,                                                    |
      | d σ_33 d F_00,                                                    |
      | ...                                                               |
      | ...                                                               |
      \ ...                                                               /

      such that it can be interpreted as a column major 6x3x3 tensor (4th order, left voigt tensor)
  */

public:
  using MarmotMaterial::MarmotMaterial;

  virtual void computeStress( double*       stress,
                              double&       K_local,
                              double&       nonLocalRadius,
                              double*       dStressDDFNew,
                              double*       dK_localDDFNew,
                              double*       dStressDK,
                              const double* FOld,
                              const double* FNew,
                              const double  KOld,
                              const double  dK,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

  virtual void computePlaneStress( double*       stress,
                                   double&       K_local,
                                   double&       nonLocalRadius,
                                   double*       dStressDDFNew,
                                   double*       dK_localDDFNew,
                                   double*       dStressDK,
                                   const double* FOld,
                                   const double* FNew,
                                   double        KOld,
                                   double        dK,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT ){};

  virtual void computeUniaxialStress( double*       stress,
                                      double&       K_local,
                                      double&       nonLocalRadius,
                                      double*       dStressDDFNew,
                                      double*       dK_localDDFNew,
                                      double*       dStressDK,
                                      const double* FOld,
                                      const double* FNew,
                                      double        KOld,
                                      double        dK,
                                      const double* timeOld,
                                      const double  dT,
                                      double&       pNewDT ){};
};
