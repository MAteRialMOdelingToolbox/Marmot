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

class MarmotMaterialHypoElastic : public MarmotMaterialMechanical {

    // Derived abstract base class for elastic materials expressed purely in rate form in terms of stretching rate d, i.e,
    // 'hypoelastic materials':
    //
    // ∇ 
    // σ = f (σ, d, t, .. ), 
    //
    // formulated incrementally as σ_np = f (σ_n, Δε, Δt, t_n, .. ) 
    // with Δε = d * Δt
    //
    // Algorithmic tangent: dσdε = d Δσ d Δε
    //
    // compatible with Abaqus interface 

  public:
    using MarmotMaterialMechanical::MarmotMaterialMechanical;

    double characteristicElementLength;
    void setCharacteristicElementLength( double length );

    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;

    // Abstract methods
    virtual void computeStress( double*       stress,
                                double*       dStressDDStrain,
                                const double* dStrain,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) = 0;

    using MarmotMaterialMechanical::computePlaneStress;
    virtual void computePlaneStress( double*       stress,
                                     double*       dStressDDStrain,
                                     double*       dStrain,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );

    using MarmotMaterialMechanical::computeUniaxialStress;
    virtual void computeUniaxialStress( double*       stress,
                                        double*       dStressDDStrain,
                                        double*       dStrain,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );
};
