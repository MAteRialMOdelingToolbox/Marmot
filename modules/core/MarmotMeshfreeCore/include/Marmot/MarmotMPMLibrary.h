/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of MaterialPoints and Structural Analysis
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
#include "Marmot/MarmotCell.h"
#include "Marmot/MarmotCellElement.h"
#include "Marmot/MarmotMaterialPoint.h"
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  // MaterialPointFactory
  //
  // - Allows materials to register themselve with their name and ID
  // - Allows the user to create instances of materials

  class MarmotMaterialPointFactory {
  public:
    using materialPointFactoryFunction = MarmotMaterialPoint* (*)( int           materialPointNumber,
                                                                   const double* vertexCoordinates,
                                                                   int           sizeVertexCoordinates,
                                                                   double        volume
                                                                   /* const MarmotMaterialSection& material */
    );
    MarmotMaterialPointFactory()       = delete;

    static MarmotMaterialPoint* createMaterialPoint( const std::string& materialPointName,
                                                     int                materialNumber,
                                                     const double*      vertexCoordinates,
                                                     int                sizeVertexCoordinates,
                                                     double             volume
                                                     /* const MarmotMaterialSection& material */
    );

    static bool registerMaterialPoint( const std::string& materialName, materialPointFactoryFunction factoryFunction );

  private:
    bool checkIfMaterialPointIsRegistered( const std::string& materialPointName );

    static std::unordered_map< std::string, materialPointFactoryFunction > materialPointFactoryFunctionByName;
  };

  // CellFactory
  //
  // - Allows cells to register themselve with their name
  // - Allows the user to create instances of cells

  class MarmotCellFactory {
  public:
    using cellFactoryFunction = MarmotCell* (*)( int           cellNumber,
                                                 const double* vertexCoordinates,
                                                 int           sizeVertexCoordinates );

    using bSplineCellFactoryFunction = MarmotCell* (*)( int           cellNumber,
                                                        const double* vertexCoordinates,
                                                        int           sizeVertexCoordinates,
                                                        const double* knotVectors,
                                                        int           sizeKnotVectors );

    MarmotCellFactory() = delete;

    static MarmotCell* createCell( const std::string& cellName,
                                   int                cellNumber,
                                   const double*      nodeCoordinates,
                                   int                sizeNodeCoordinates );

    static MarmotCell* createBSplineCell( const std::string& cellName,
                                          int                cellNumber,
                                          const double*      nodeCoordinates,
                                          int                sizeNodeCoordinates,
                                          const double*      knotVectors,
                                          int                sizeKnotVectors );

    static bool registerCell( const std::string& cellName, cellFactoryFunction factoryFunction );
    static bool registerBSplineCell( const std::string& cellName, bSplineCellFactoryFunction factoryFunction );

  private:
    static std::unordered_map< std::string, cellFactoryFunction >        cellFactoryFunctionByName;
    static std::unordered_map< std::string, bSplineCellFactoryFunction > bSplineCellFactoryFunctionByName;
  };

  // CellElementFactory
  //
  // - Allows cellElements to register themselve with their name
  // - Allows the user to create instances of cellElements

  class MarmotCellElementFactory {
  public:
    using cellElementFactoryFunction = MarmotCellElement* (*)( int                cellElementNumber,
                                                               const double*      vertexCoordinates,
                                                               int                sizeVertexCoordinates,
                                                               const std::string& quadratureRule,
                                                               int                quadratureOrder );

    MarmotCellElementFactory() = delete;

    static MarmotCellElement* createCellElement( const std::string& cellElementName,
                                                 int                cellElementNumber,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 const std::string& quadratureRule,
                                                 int                quadratureOrder );

    static bool registerCellElement( const std::string& cellElementName, cellElementFactoryFunction factoryFunction );

    /* private: */
    static std::unordered_map< std::string, cellElementFactoryFunction > cellElementFactoryFunctionByName;
  };

} // namespace MarmotLibrary
