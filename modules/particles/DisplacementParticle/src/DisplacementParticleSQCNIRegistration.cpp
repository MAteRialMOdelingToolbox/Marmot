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
#include "Marmot/DisplacementParticleSQCNI.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace {
  using namespace Marmot::Discretizations::Meshfree;


  const static bool DisplacementParticleSQCNI_PlaneStrain_Quad_isRegistered = Marmot::Registry::ParticleFactory::
    registerParticle( "DisplacementSQCNI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Discretizations::Meshfree::Approximation& approximation )
                        -> Marmot::Discretizations::Meshfree::Particle* {
                        return new DisplacementParticleSQCNI<
                          2,
                          4 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNI< 2, 4 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementParticleSQCNI_3D_Hexa_isRegistered = Marmot::Registry::ParticleFactory::
    registerParticle( "DisplacementSQCNI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Discretizations::Meshfree::Approximation& approximation )
                        -> Marmot::Discretizations::Meshfree::Particle* {
                        return new DisplacementParticleSQCNI<
                          3,
                          8 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNI< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementParticleSNNI_PlaneStrain_Quad_isRegistered = Marmot::Registry::ParticleFactory::
    registerParticle( "DisplacementSNNI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Discretizations::Meshfree::Approximation& approximation )
                        -> Marmot::Discretizations::Meshfree::Particle* {
                        return new DisplacementParticleSQCNI< 2, 4 >( cellID,
                                                                      nodeCoordinates,
                                                                      sizeNodeCoordinates,
                                                                      volume,
                                                                      materialName,
                                                                      materialProperties,
                                                                      sizeMaterialProperties,
                                                                      approximation,
                                                                      DisplacementParticleSQCNI< 2, 4 >::
                                                                        SmoothingDomainUpdateType::None );
                      } );

  const static bool DisplacementParticleSQCNI_R_PlaneStrain_Quad_isRegistered = Marmot::Registry::ParticleFactory::
    registerParticle( "DisplacementSQCNI_R/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Discretizations::Meshfree::Approximation& approximation )
                        -> Marmot::Discretizations::Meshfree::Particle* {
                        return new DisplacementParticleSQCNI< 2, 4 >( cellID,
                                                                      nodeCoordinates,
                                                                      sizeNodeCoordinates,
                                                                      volume,
                                                                      materialName,
                                                                      materialProperties,
                                                                      sizeMaterialProperties,
                                                                      approximation,
                                                                      DisplacementParticleSQCNI< 2, 4 >::
                                                                        SmoothingDomainUpdateType::RotationOnly );
                      } );

  const static bool DisplacementParticleSQCNI_RU_PlaneStrain_Quad_isRegistered = Marmot::Registry::ParticleFactory::
    registerParticle( "DisplacementSQCNI_RU/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Discretizations::Meshfree::Approximation& approximation )
                        -> Marmot::Discretizations::Meshfree::Particle* {
                        return new DisplacementParticleSQCNI<
                          2,
                          4 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNI< 2,
                                                          4 >::SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                      } );

} // namespace
