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
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotParticle.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>

namespace Marmot::Registry {

  std::string makeStringUpperCase_( const std::string& str )
  {
    std::string strUpperCase = str;

    std::transform( strUpperCase.begin(), strUpperCase.end(), strUpperCase.begin(), toupper );

    return strUpperCase;
  }

  std::unordered_map< std::string, ParticleFactory::particleFactoryFunction >& ParticleFactory::
    particleFactoryFunctionByName()
  {
    static std::unordered_map< std::string, particleFactoryFunction > map;
    return map;
  }

  bool ParticleFactory::registerParticle( const std::string& particleName, particleFactoryFunction factoryFunction )
  {

    const auto particleNameUpperCase = makeStringUpperCase_( particleName );

    assert( particleFactoryFunctionByName().find( particleNameUpperCase ) == particleFactoryFunctionByName().end() );

    particleFactoryFunctionByName()[particleNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Discretizations::Meshfree::Particle* ParticleFactory::createParticle(
    const std::string& particleName,
    int                particleNumber,
    const double*      vertexCoordinates,
    int                sizeVertexCoordinates,
    double             volume,
    // Marmot::Discretizations::MaterialPoints::MaterialPoint&                                 mp,
    const std::string&                                      materialName,
    const double*                                           materialProperties,
    int                                                     sizeMaterialProperties,
    const Marmot::Discretizations::Meshfree::Approximation& approximation )
  {
    const auto particleNameUpperCase = makeStringUpperCase_( particleName );

    try {
      return particleFactoryFunctionByName().at( particleNameUpperCase )( particleNumber,
                                                                          vertexCoordinates,
                                                                          sizeVertexCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << "Invalid particle " << particleName << " requested!" );
    }
  }

} // namespace Marmot::Registry
