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
#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace {
  using namespace Marmot::Discretizations::MaterialPoints;

  /**
   * @brief Register a gradient-enhanced finite-strain material point.
   * @tparam nDim      Spatial dimension.
   * @tparam nNonlocal Number of nonlocal variables.
   * @param name The registered name.
   * @return True if the registration succeeded.
   */
  template < int nDim, int nNonlocal >
  bool registerGradientEnhancedFiniteStrainMaterialPoint( const std::string& name )
  {
    return Marmot::Registry::MaterialPointFactory::
      registerMaterialPoint( name,
                             []( int mpNumber, const double* vertexCoordinates, int nVertexCoordinates, double volume )
                               -> Marmot::Discretizations::MaterialPoints::MaterialPoint* {
                               return new GradientEnhancedFiniteStrainMaterialPoint< nDim,
                                                                                     nNonlocal >( mpNumber,
                                                                                                  vertexCoordinates,
                                                                                                  nVertexCoordinates,
                                                                                                  volume );
                             } );
  }

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint2D_isRegistered = registerGradientEnhancedFiniteStrainMaterialPoint< 2,
                                                                                                                  1 >(
      "GradientEnhancedFiniteStrain/PlaneStrain" );

  const static bool
    GradientEnhancedFiniteStrainMaterialPoint3D_isRegistered = registerGradientEnhancedFiniteStrainMaterialPoint< 3,
                                                                                                                  1 >(
      "GradientEnhancedFiniteStrain/3D" );

  const static bool
    GradientEnhanced6FiniteStrainMaterialPoint3D_isRegistered = registerGradientEnhancedFiniteStrainMaterialPoint< 3,
                                                                                                                   6 >(
      "GradientEnhanced6FiniteStrain/3D" );

} // namespace
