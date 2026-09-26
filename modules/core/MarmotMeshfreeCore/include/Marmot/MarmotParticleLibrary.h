/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Particles and Structural Analysis
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
// #include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotParticle.h"
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  // ParticleFactory
  //
  // - Allows materials to register themselve with their name and ID
  // - Allows the user to create instances of materials

  class MarmotParticleFactory {
  public:
    using particleFactoryFunction =
      Marmot::Meshfree::MarmotParticle* (*)( int           particleNumber,
                                             const double* vertexCoordinates,
                                             int           sizeVertexCoordinates,
                                             double        volume,
                                             // MarmotMaterialPoint&                                 mp,
                                             //
                                             const std::string& materialName,
                                             const double*      materialProperties,
                                             int                nMaterialProperties,

                                             const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation );
    MarmotParticleFactory() = delete;

    static Marmot::Meshfree::MarmotParticle* createParticle(
      const std::string& particleName,
      int                materialNumber,
      const double*      vertexCoordinates,
      int                sizeVertexCoordinates,
      double             volume,
      // MarmotMaterialPoint&                                 mp,
      //
      const std::string& materialName,
      const double*      materialProperties,
      int                nMaterialProperties,

      const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation );

    static bool registerParticle( const std::string& particleName, particleFactoryFunction factoryFunction );

  private:
    bool checkIfParticleIsRegistered( const std::string& particleName );

    static std::unordered_map< std::string, particleFactoryFunction > particleFactoryFunctionByName;
  };

} // namespace MarmotLibrary
