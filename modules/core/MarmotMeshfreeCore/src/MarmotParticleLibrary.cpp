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
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotParticle.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  std::string makeStringUpperCase_( const std::string& str )
  {
    std::string strUpperCase = str;

    std::transform( strUpperCase.begin(), strUpperCase.end(), strUpperCase.begin(), toupper );

    return strUpperCase;
  }

  std::unordered_map< std::string, MarmotParticleFactory::particleFactoryFunction >
    MarmotParticleFactory::particleFactoryFunctionByName;

  bool MarmotParticleFactory::registerParticle( const std::string&      particleName,
                                                particleFactoryFunction factoryFunction )
  {

    const auto particleNameUpperCase = makeStringUpperCase_( particleName );

    assert( particleFactoryFunctionByName.find( particleNameUpperCase ) == particleFactoryFunctionByName.end() );

    particleFactoryFunctionByName[particleNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Meshfree::MarmotParticle* MarmotParticleFactory::createParticle(
    const std::string& particleName,
    int                particleNumber,
    const double*      vertexCoordinates,
    int                sizeVertexCoordinates,
    double             volume,
    // MarmotMaterialPoint&                                 mp,
    const std::string&                                   materialName,
    const double*                                        materialProperties,
    int                                                  sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
  {
    const auto particleNameUpperCase = makeStringUpperCase_( particleName );

    try {
      return particleFactoryFunctionByName.at( particleNameUpperCase )( particleNumber,
                                                                        vertexCoordinates,
                                                                        sizeVertexCoordinates,
                                                                        volume,
                                                                        materialName,
                                                                        materialProperties,
                                                                        sizeMaterialProperties,
                                                                        approximation );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid particle " << particleName << " requested!" );
    }
  }

} // namespace MarmotLibrary
