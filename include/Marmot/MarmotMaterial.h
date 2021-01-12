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
#include "Marmot/MarmotUtils.h"
#include <string>

class MarmotMaterial {

  protected:

    const double* materialProperties;
    const int nMaterialProperties;

    double* stateVars;
    int     nStateVars;

  public:
    const int materialNumber;

    MarmotMaterial( const double*materialProperties, int nMaterialProperties, int materialNumber );

    virtual ~MarmotMaterial();

    virtual int getNumberOfRequiredStateVars() = 0;

    virtual void assignStateVars( double* stateVars, int nStateVars );

    virtual PermanentResultLocation getPermanentResultPointer( const std::string& resultName ) = 0;

    const double* getAssignedStateVars();

    int getNumberOfAssignedStateVars();
};
