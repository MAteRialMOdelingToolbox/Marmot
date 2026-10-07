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
#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace Marmot::MaterialPoints::Registration {

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint2D_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
      registerMaterialPoint( "GradientEnhancedFiniteStrain/PlaneStrain",
                             []( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
                               -> MarmotMaterialPoint* {
                               return new GradientEnhancedFiniteStrainMaterialPoint2D( mpNumber,
                                                                                       vertexCoordinates,
                                                                                       nVertexCoordinates,
                                                                                       volume );
                             } );

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint3D_isRegistered = MarmotLibrary::MarmotMaterialPointFactory::
      registerMaterialPoint( "GradientEnhancedFiniteStrain/3D",
                             []( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
                               -> MarmotMaterialPoint* {
                               return new GradientEnhancedFiniteStrainMaterialPoint3D( mpNumber,
                                                                                       vertexCoordinates,
                                                                                       nVertexCoordinates,
                                                                                       volume );
                             } );

} // namespace Marmot::MaterialPoints::Registration
