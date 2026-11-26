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
#include "Marmot/GradientEnhancedFiniteStrainParticle.h"

namespace {
  using namespace Marmot::Discretizations::Meshfree;

  const static bool
    GradientEnhancedFiniteStrainParticle_PlaneStrain_isRegistered = registerGradientEnhancedFiniteStrainParticle<
      GradientEnhancedFiniteStrainParticle< 2, 1 > >( "GradientEnhancedFiniteStrain/PlaneStrain/Point" );

  const static bool GradientEnhancedFiniteStrainParticle_3D_isRegistered = registerGradientEnhancedFiniteStrainParticle<
    GradientEnhancedFiniteStrainParticle< 3, 1 > >( "GradientEnhancedFiniteStrain/3D/Point" );

  const static bool
    GradientEnhanced6FiniteStrainParticle_3D_isRegistered = registerGradientEnhancedFiniteStrainParticle<
      GradientEnhancedFiniteStrainParticle< 3, 6 > >( "GradientEnhanced6FiniteStrain/3D/Point" );

} // namespace
