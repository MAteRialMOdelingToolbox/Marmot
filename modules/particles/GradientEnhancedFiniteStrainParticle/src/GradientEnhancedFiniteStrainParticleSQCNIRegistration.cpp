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
#include "Marmot/GradientEnhancedFiniteStrainParticleRegistration.h"
#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNI.h"

namespace {
  using namespace Marmot::Discretizations::Meshfree;

  const static bool
    GradientEnhancedFiniteStrainSQCNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSNNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >::SmoothingDomainUpdateType::None >(
      "GradientEnhancedFiniteStrainSNNI/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_R_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >::SmoothingDomainUpdateType::RotationOnly >(
      "GradientEnhancedFiniteStrainSQCNI_R/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RU_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 2, 4, 1 >::SmoothingDomainUpdateType::RotationAndPrincipalStretch >(
      "GradientEnhancedFiniteStrainSQCNI_RU/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhancedFiniteStrainSQCNI/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSNNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >::SmoothingDomainUpdateType::None >(
      "GradientEnhancedFiniteStrainSNNI/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_R_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >::SmoothingDomainUpdateType::RotationOnly >(
      "GradientEnhancedFiniteStrainSQCNI_R/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RU_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 1 >::SmoothingDomainUpdateType::RotationAndPrincipalStretch >(
      "GradientEnhancedFiniteStrainSQCNI_RU/3D/Hexa" );

  const static bool
    GradientEnhanced6FiniteStrainSQCNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 6 >,
      GradientEnhancedFiniteStrainParticleSQCNI< 3, 8, 6 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhanced6FiniteStrainSQCNI/3D/Hexa" );

} // namespace
