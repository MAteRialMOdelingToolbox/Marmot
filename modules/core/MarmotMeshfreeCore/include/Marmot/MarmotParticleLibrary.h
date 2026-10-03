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

  /**
   * @class MarmotLibrary::MarmotParticleFactory
   * @brief Registry and factory of the meshfree particles (Marmot::Meshfree::MarmotParticle).
   *
   * - Allows particles to register themselves with their name (see registerParticle(); particle modules do this in
   *   their *Registration.cpp files at static initialization).
   * - Allows the user (the host framework) to create instances of particles by name (createParticle()).
   *
   * Particle names are case-insensitive: they are converted to upper case both at registration and at creation.
   * They follow the pattern `<Physics><IntegrationScheme>/<Dimension>/<Shape>`, e.g.,
   * @c DisplacementSQCNI/PlaneStrain/Quad or @c GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa.
   */
  class MarmotParticleFactory {
  public:
    /**
     * @brief Signature of the factory function of a particle type.
     * @details Arguments: particle number, vertex coordinates (the center for point particles, the cell vertices
     * otherwise) and their number, the volume (for cell-shaped particles zero, the volume follows from the vertices),
     * the material name, the material properties and their number, and the meshfree approximation.
     */
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
    /// @brief The factory is a static class and cannot be instantiated.
    MarmotParticleFactory() = delete;

    /**
     * @brief Create a particle by its registered name.
     * @param[in] particleName Registered name of the particle type (case-insensitive).
     * @param[in] materialNumber The number (label) of the particle, passed as the first argument of the factory
     * function.
     * @param[in] vertexCoordinates Vertex coordinates (the center for point particles), vertex by vertex.
     * @param[in] sizeVertexCoordinates Number of values in @p vertexCoordinates.
     * @param[in] volume Volume of the particle (point particles), or zero for cell-shaped particles.
     * @param[in] materialName Name of the material.
     * @param[in] materialProperties Material properties.
     * @param[in] nMaterialProperties Number of material properties.
     * @param[in] approximation The meshfree approximation; it must outlive the particle.
     * @return A new particle, owned by the caller.
     * @throws std::invalid_argument if no particle is registered under @p particleName.
     */
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

    /**
     * @brief Register a particle type under a name.
     * @details A name must be registered only once (checked by an assertion in debug builds).
     * @param[in] particleName Name of the particle type (stored in upper case).
     * @param[in] factoryFunction Function creating an instance.
     * @return true, to allow the registration in the initializer of a static variable.
     */
    static bool registerParticle( const std::string& particleName, particleFactoryFunction factoryFunction );

  private:
    /**
     * @brief Check whether a particle name is registered (declared only, not defined).
     * @param[in] particleName Name of the particle type.
     * @return true if registered.
     */
    bool checkIfParticleIsRegistered( const std::string& particleName );

    /// @brief Registered factory functions by upper-case particle name (a function-local static: registrations run
    /// during static initialization, possibly before a static data member of this translation unit would be
    /// constructed).
    static std::unordered_map< std::string, particleFactoryFunction >& particleFactoryFunctionByName();
  };

} // namespace MarmotLibrary
