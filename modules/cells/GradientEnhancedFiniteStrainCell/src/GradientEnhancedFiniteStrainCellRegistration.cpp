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
#include "Marmot/GradientEnhancedFiniteStrainCell.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace Marmot::Cells::Registration {

  using namespace MarmotLibrary;

  const static bool CPE4CELL_isRegistered = MarmotLibrary::MarmotCellFactory::
    registerCell( "GradientEnhancedFiniteStrain/Quad4",
                  []( int cellID, const double* nodeCoordinates, int sizeNodeCoordinates ) -> MarmotCell* {
                    return new LagrangianGradientEnhancedFiniteStrainCell< 2, 4 >( cellID,
                                                                                   nodeCoordinates,
                                                                                   sizeNodeCoordinates );
                  } );

  const static bool C3D8CELL_isRegistered = MarmotLibrary::MarmotCellFactory::
    registerCell( "GradientEnhancedFiniteStrain/Hexa8",
                  []( int cellID, const double* nodeCoordinates, int sizeNodeCoordinates ) -> MarmotCell* {
                    return new LagrangianGradientEnhancedFiniteStrainCell< 3, 8 >( cellID,
                                                                                   nodeCoordinates,
                                                                                   sizeNodeCoordinates );
                  } );

  const static bool BSpline2D_1 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/1",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 2, 4, 1 >( cellID,
                                                                                          nodeCoordinates,
                                                                                          sizeNodeCoordinates,
                                                                                          knotVectors,
                                                                                          sizeKnotVectors );
                         } );

  const static bool BSpline2D_2 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/2",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 2, 9, 2 >( cellID,
                                                                                          nodeCoordinates,
                                                                                          sizeNodeCoordinates,
                                                                                          knotVectors,
                                                                                          sizeKnotVectors );
                         } );

  const static bool BSpline2D_3 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/3",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 2, 16, 3 >( cellID,
                                                                                           nodeCoordinates,
                                                                                           sizeNodeCoordinates,
                                                                                           knotVectors,
                                                                                           sizeKnotVectors );
                         } );

  const static bool BSpline3D_1 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/3D/1",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 3, 2 * 2 * 2, 1 >( cellID,
                                                                                                  nodeCoordinates,
                                                                                                  sizeNodeCoordinates,
                                                                                                  knotVectors,
                                                                                                  sizeKnotVectors );
                         } );

  const static bool BSpline3D_2 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/3D/2",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 3, 3 * 3 * 3, 2 >( cellID,
                                                                                                  nodeCoordinates,
                                                                                                  sizeNodeCoordinates,
                                                                                                  knotVectors,
                                                                                                  sizeKnotVectors );
                         } );

  const static bool BSpline3D_3 = MarmotLibrary::MarmotCellFactory::
    registerBSplineCell( "GradientEnhancedFiniteStrain/BSpline/3D/3",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> MarmotCell* {
                           return new BSplineGradientEnhancedFiniteStrainCell< 3, 4 * 4 * 4, 3 >( cellID,
                                                                                                  nodeCoordinates,
                                                                                                  sizeNodeCoordinates,
                                                                                                  knotVectors,
                                                                                                  sizeKnotVectors );
                         } );

} // namespace Marmot::Cells::Registration
