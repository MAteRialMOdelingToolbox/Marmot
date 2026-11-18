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
#include "Marmot/DisplacementFBarCell.h"
#include "Marmot/MarmotMPMLibrary.h"

namespace {
  using namespace Marmot::Discretizations::Cells;

  const static bool CPE4FBARCELL_isRegistered = Marmot::Registry::CellFactory::
    registerCell( "DisplacementFBar/Quad4",
                  []( int           cellID,
                      const double* nodeCoordinates,
                      int           sizeNodeCoordinates ) -> Marmot::Discretizations::Cells::Cell* {
                    return new LagrangianDisplacementFBarCell< 2, 4 >( cellID, nodeCoordinates, sizeNodeCoordinates );
                  } );

  const static bool C3D8FBARCELL_isRegistered = Marmot::Registry::CellFactory::
    registerCell( "DisplacementFBar/Hexa8",
                  []( int           cellID,
                      const double* nodeCoordinates,
                      int           sizeNodeCoordinates ) -> Marmot::Discretizations::Cells::Cell* {
                    return new LagrangianDisplacementFBarCell< 3, 8 >( cellID, nodeCoordinates, sizeNodeCoordinates );
                  } );

  const static bool BSpline2D_1 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/1",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 2, 4, 1 >( cellID,
                                                                              nodeCoordinates,
                                                                              sizeNodeCoordinates,
                                                                              knotVectors,
                                                                              sizeKnotVectors );
                         } );

  const static bool BSpline2D_2 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/2",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 2, 9, 2 >( cellID,
                                                                              nodeCoordinates,
                                                                              sizeNodeCoordinates,
                                                                              knotVectors,
                                                                              sizeKnotVectors );
                         } );

  const static bool BSpline2D_3 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/3",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 2, 16, 3 >( cellID,
                                                                               nodeCoordinates,
                                                                               sizeNodeCoordinates,
                                                                               knotVectors,
                                                                               sizeKnotVectors );
                         } );

  const static bool BSpline3D_1 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/3D/1",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 3, 2 * 2 * 2, 1 >( cellID,
                                                                                      nodeCoordinates,
                                                                                      sizeNodeCoordinates,
                                                                                      knotVectors,
                                                                                      sizeKnotVectors );
                         } );

  const static bool BSpline3D_2 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/3D/2",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 3, 3 * 3 * 3, 2 >( cellID,
                                                                                      nodeCoordinates,
                                                                                      sizeNodeCoordinates,
                                                                                      knotVectors,
                                                                                      sizeKnotVectors );
                         } );

  const static bool BSpline3D_3 = Marmot::Registry::CellFactory::
    registerBSplineCell( "DisplacementFBar/BSpline/3D/3",
                         []( int           cellID,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors ) -> Marmot::Discretizations::Cells::Cell* {
                           return new BSplineDisplacementFBarCell< 3, 4 * 4 * 4, 3 >( cellID,
                                                                                      nodeCoordinates,
                                                                                      sizeNodeCoordinates,
                                                                                      knotVectors,
                                                                                      sizeKnotVectors );
                         } );

} // namespace
