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

class MarmotMaterialHyperElastic : public MarmotMaterialMechanical {

    // Derived abstract base class for a _simple_, purely hyperelastic material to be used within TL elements
    //
    // stress measure: Piola - Kirchhoff II .. S
    // strain measure for algorithmic tangent: Green - Lagrange .. E = 1/2 ( F^T * F - I )
    //
    //      ∂ f( E, t )
    // S =  -----------
    //      ∂    E 
    //
    // Algorithmic tangent: dS/dE

  public:
    using MarmotMaterialMechanical::MarmotMaterialMechanical;

    // Default implementation provided
    virtual void computeStress( double*       S,    // PK2
                                double*       dSdE, // d PK2 d GL_E
                                const double* FOld,
                                const double* FNew,
                                const double* timeOld,
                                const double  dT,
                                double&       pNewDT ) override;

    // Abstract methods
    virtual void computeStressPK2( double*       S,    // PK2
                                   double*       dSdE, // d PK2 d GL_E
                                   const double* E,
                                   const double* timeOld,
                                   const double  dT,
                                   double&       pNewDT ) = 0;

    virtual void computePlaneStressPK2( double*       S,
                                        double*       dSdE,
                                        double*       E,
                                        const double* timeOld,
                                        const double  dT,
                                        double&       pNewDT );

    virtual void computeUniaxialStressPK2( double*       S,
                                           double*       dSdE,
                                           double*       E,
                                           const double* timeOld,
                                           const double  dT,
                                           double&       pNewDT );
};
