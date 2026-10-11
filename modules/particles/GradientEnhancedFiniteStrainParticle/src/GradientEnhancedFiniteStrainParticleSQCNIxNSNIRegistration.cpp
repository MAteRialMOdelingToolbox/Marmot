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
#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNIxNSNI.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool
    GradientEnhancedFiniteStrainParticleSNNIxNSNI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrainSNNIxNSNI/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2,
                                                                                 4 >::SmoothingDomainUpdateType::None );
                        } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_R_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_RU_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::DeformationGradient );
                                             } );

  const static bool
    GradientEnhancedFiniteStrainParticleSNNIxNSNI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "GradientEnhancedFiniteStrainSNNIxNSNI/3D/Hexa",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                            3,
                            8 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3,
                                                                                 8 >::SmoothingDomainUpdateType::None );
                        } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationOnly );
                                             } );

  const static bool GradientEnhancedFiniteStrainParticleSQCNI_RUxNSNI_3D_Hexa_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/3D/Hexa",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new GradientEnhancedFiniteStrainParticleSQCNIxNSNI<
                                                 3,
                                                 8 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );

} // namespace Marmot::Meshfree
