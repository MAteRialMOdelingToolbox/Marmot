/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
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
#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace Marmot::MaterialPoints::Registration {

  const static bool DisplacementPlaneStrainMaterialPoint_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
    registerMaterialPoint( "Displacement/PlaneStrain",
                           []( int           materialPointNumber,
                               const double* vertexCoordinates,
                               int           nVertexCoordinates,
                               double        volume ) -> MarmotMaterialPoint* {
                             return new DisplacementMaterialPoint2D( materialPointNumber,
                                                                     vertexCoordinates,
                                                                     nVertexCoordinates,
                                                                     volume );
                           } );

  const static bool Displacement3DMaterialPoint_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
    registerMaterialPoint( "Displacement/3D",
                           []( int           materialPointNumber,
                               const double* vertexCoordinates,
                               int           nVertexCoordinates,
                               double        volume ) -> MarmotMaterialPoint* {
                             return new DisplacementMaterialPoint3D( materialPointNumber,
                                                                     vertexCoordinates,
                                                                     nVertexCoordinates,
                                                                     volume );
                           } );

} // namespace Marmot::MaterialPoints::Registration
