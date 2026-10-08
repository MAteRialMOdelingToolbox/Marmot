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
#include "Marmot/DisplacementParticleSQCNIxSDI.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool DisplacementParticleSQCNIxSDI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementSQCNIxSDI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxSDI<
                          2,
                          4 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNIxSDI< 2, 4 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementParticleSQCNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementSQCNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxSDI<
                          3,
                          8 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementParticleSQCNIxSDI< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementParticleRS_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementRS-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::
                                                                              RotationAndPrincipalStretch );
                      } );

  const static bool DisplacementParticleR_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementR-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::RotationOnly );
                      } );

  const static bool DisplacementParticleSNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementSNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::None );
                      } );

} // namespace Marmot::Meshfree
