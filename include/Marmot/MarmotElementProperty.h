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

namespace MarmotLibrary {
    enum MaterialCode : int;
}

class MarmotMaterialSection  {
  public:
    MarmotLibrary::MaterialCode materialCode;
    const double*             materialProperties;
    int                       nMaterialProperties;

    MarmotMaterialSection( MarmotLibrary::MaterialCode materialCode,
                        const double*             materialProperties,
                        int                       nMaterialProperties )
        : materialCode( materialCode ),
          materialProperties( materialProperties ),
          nMaterialProperties( nMaterialProperties ){};
};

class ElementProperties  {
  public:
    const double* elementProperties;
    int           nElementProperties;

    ElementProperties( const double* elementProperties, int nElementProperties )
        : elementProperties( elementProperties ),
          nElementProperties( nElementProperties ){};
};
