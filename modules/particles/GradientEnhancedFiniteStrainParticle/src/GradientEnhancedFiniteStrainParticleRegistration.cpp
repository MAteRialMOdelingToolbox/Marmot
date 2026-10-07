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
#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool
    GradientEnhancedFiniteStrainParticle_PlaneStrain_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrain/PlaneStrain/Point",
                        []( int           cellID,
                            const double* nodeCoordinates,
                            int           sizeNodeCoordinates,
                            double        volume,
                            // MarmotMaterialPoint&                                 mp,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticle< 2 >( cellID,
                                                                                nodeCoordinates,
                                                                                sizeNodeCoordinates,
                                                                                volume,
                                                                                // mp,
                                                                                materialName,
                                                                                materialProperties,
                                                                                sizeMaterialProperties,
                                                                                approximation );
                        } );

  const static bool GradientEnhancedFiniteStrainParticle_3D_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "GradientEnhancedFiniteStrain/3D/Point",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new GradientEnhancedFiniteStrainParticle< 3 >( cellID,
                                                                              nodeCoordinates,
                                                                              sizeNodeCoordinates,
                                                                              volume,
                                                                              materialName,
                                                                              materialProperties,
                                                                              sizeMaterialProperties,
                                                                              approximation );
                      } );

} // namespace Marmot::Meshfree
