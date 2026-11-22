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
#pragma once

#include "Marmot/MarmotMeshfreeApproximation.h"
#include "Marmot/MarmotParticle.h"
#include "Marmot/MarmotParticleLibrary.h"
#include <string>

namespace Marmot::Discretizations::Meshfree {

  /**
   * @brief Register a gradient-enhanced finite-strain point particle with the Marmot::Registry::ParticleFactory.
   * @tparam ParticleType A GradientEnhancedFiniteStrainParticle.
   * @param[in] name The registered name.
   * @return True if the registration succeeded.
   */
  template < class ParticleType >
  bool registerGradientEnhancedFiniteStrainParticle( const std::string& name )
  {
    return Marmot::Registry::ParticleFactory::registerParticle( name,
                                                                []( int                  particleID,
                                                                    const double*        coordinates,
                                                                    int                  sizeCoordinates,
                                                                    double               volume,
                                                                    const std::string&   materialName,
                                                                    const double*        materialProperties,
                                                                    int                  sizeMaterialProperties,
                                                                    const Approximation& approximation ) -> Particle* {
                                                                  return new ParticleType( particleID,
                                                                                           coordinates,
                                                                                           sizeCoordinates,
                                                                                           volume,
                                                                                           materialName,
                                                                                           materialProperties,
                                                                                           sizeMaterialProperties,
                                                                                           approximation );
                                                                } );
  }

  /**
   * @brief Register a gradient-enhanced finite-strain particle with a smoothing domain (SQCNI, SNNI, ... and their
   * NSNI variants) with the Marmot::Registry::ParticleFactory.
   * @tparam ParticleType A GradientEnhancedFiniteStrainParticleSQCNI or a class derived from it.
   * @tparam updateType   How the smoothing domain follows the deformation.
   * @param[in] name The registered name.
   * @return True if the registration succeeded.
   */
  template < class ParticleType, typename ParticleType::SmoothingDomainUpdateType updateType >
  bool registerGradientEnhancedFiniteStrainParticleWithDomain( const std::string& name )
  {
    return Marmot::Registry::ParticleFactory::registerParticle( name,
                                                                []( int                  particleID,
                                                                    const double*        vertexCoordinates,
                                                                    int                  sizeVertexCoordinates,
                                                                    double               volume,
                                                                    const std::string&   materialName,
                                                                    const double*        materialProperties,
                                                                    int                  sizeMaterialProperties,
                                                                    const Approximation& approximation ) -> Particle* {
                                                                  return new ParticleType( particleID,
                                                                                           vertexCoordinates,
                                                                                           sizeVertexCoordinates,
                                                                                           volume,
                                                                                           materialName,
                                                                                           materialProperties,
                                                                                           sizeMaterialProperties,
                                                                                           approximation,
                                                                                           updateType );
                                                                } );
  }

} // namespace Marmot::Discretizations::Meshfree
