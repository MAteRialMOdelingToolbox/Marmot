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

class MarmotMaterialMechanical : public MarmotMaterial {

  /*
     Abstract basic class for Mechanical materials with scalar nonlocal interaction.

     Formulated incrementally as σ_np, K_local_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, , Kn, ΔK, K_local_n.. )

     Algorithmic tangents: dσdF = d σ_np d (dxdX_np)
                           dK_LocaldF = d K_local_np d (dxdX_np)
                           dσdK = d σ_np d ΔK
  */

public:
  using MarmotMaterial::MarmotMaterial;

  virtual void computeStress( double*       stress,
                              double*       dStressDDFNew,
                              const double* FOld,
                              const double* FNew,
                              const double* timeOld,
                              const double  dT,
                              double&       pNewDT ) = 0;

  virtual void computePlaneStress( double*       stress,
                                   double*       dStressDDFNew,
                                   const double* FOld,
                                   double*       FNew,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT );

  virtual void computeUniaxialStress( double*       stress,
                                      double*       dStressDDFNew,
                                      const double* FOld,
                                      double*       FNew,
                                      const double* timeOld,
                                      const double  dT,
                                      double&       pNewDT );
};
