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
#include "Marmot/GradientEnhancedFiniteStrainParticleSQCNIxNSNI.h"

namespace {
  using namespace Marmot::Discretizations::Meshfree;

  const static bool
    GradientEnhancedFiniteStrainSQCNIxNSNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSNNIxNSNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >::SmoothingDomainUpdateType::None >(
      "GradientEnhancedFiniteStrainSNNIxNSNI/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RxNSNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >::SmoothingDomainUpdateType::RotationOnly >(
      "GradientEnhancedFiniteStrainSQCNI_RxNSNI/PlaneStrain/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RUxNSNI_PlaneStrain_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1 >::SmoothingDomainUpdateType::
        RotationAndPrincipalStretch >( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/PlaneStrain/Quad" );

  // the axisymmetric twin of the plane-strain particle, in the meridian half plane (r, z); with the update
  // DeformationGradient only, whose smoothing domains tile the body (see ParticleSectionIntegrals)
  using AxisymmetricSQCNIxNSNI = GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 2, 4, 1, Section::Axisymmetric >;

  const static bool
    GradientEnhancedFiniteStrainSQCNIxNSNI_Axisymmetric_Quad_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      AxisymmetricSQCNIxNSNI,
      AxisymmetricSQCNIxNSNI::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhancedFiniteStrainSQCNIxNSNI/Axisymmetric/Quad" );

  const static bool
    GradientEnhancedFiniteStrainSQCNIxNSNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSNNIxNSNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >::SmoothingDomainUpdateType::None >(
      "GradientEnhancedFiniteStrainSNNIxNSNI/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RxNSNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >::SmoothingDomainUpdateType::RotationOnly >(
      "GradientEnhancedFiniteStrainSQCNI_RxNSNI/3D/Hexa" );

  const static bool
    GradientEnhancedFiniteStrainSQCNI_RUxNSNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 1 >::SmoothingDomainUpdateType::
        RotationAndPrincipalStretch >( "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/3D/Hexa" );

  const static bool
    GradientEnhanced6FiniteStrainSQCNIxNSNI_3D_Hexa_isRegistered = registerGradientEnhancedFiniteStrainParticleWithDomain<
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 6 >,
      GradientEnhancedFiniteStrainParticleSQCNIxNSNI< 3, 8, 6 >::SmoothingDomainUpdateType::DeformationGradient >(
      "GradientEnhanced6FiniteStrainSQCNIxNSNI/3D/Hexa" );

} // namespace
