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
#include "Marmot/MarmotMaterialGradientEnhancedMechanical.h"

class MarmotMaterialGradientEnhancedHypoElastic : public MarmotMaterialGradientEnhancedMechanical {

  public:
    using MarmotMaterialGradientEnhancedMechanical::MarmotMaterialGradientEnhancedMechanical;
    
    // Abstract methods
    virtual void computeStress( double*       stress,
                                double&       K_local,
                                double&       nonLocalRadius,
                                double*       dStressDDDeformationGradient,
                                double*       dK_localDDeformationGradient,
                                double*       dStressDK,
                                const double* FOld,
                                const double* FNew,
                                const double  KOld,
                                const double  dK,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;


    // Abstract methods
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

    using MarmotMaterialGradientEnhancedMechanical::computePlaneStress;
    virtual void computePlaneStress( double*       stress,
                                     double&       K_local,
                                     double&       nonLocalRadius,
                                     double*       dStressDDStrain,
                                     double*       dK_localDDStrain,
                                     double*       dStressDK,
                                     double*       dStrain,
                                     double        KOld,
                                     double        dK,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT );
};
