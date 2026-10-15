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
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotCell.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialPoint.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>

namespace Marmot::Registry {

  std::string makeStringUpperCase( const std::string& str )
  {
    std::string strUpperCase = str;

    // TODO: MOVE TO MarmotUtils or similar

    std::transform( strUpperCase.begin(), strUpperCase.end(), strUpperCase.begin(), toupper );

    return strUpperCase;
  }

  /*
   *  __  __       _            _       _             _       _
   * |  \/  | __ _| |_ ___ _ __(_) __ _| |_ __   ___ (_)_ __ | |_ ___
   * | |\/| |/ _` | __/ _ \ '__| |/ _` | | '_ \ / _ \| | '_ \| __/ __|
   * | |  | | (_| | ||  __/ |  | | (_| | | |_) | (_) | | | | | |_\__ \
   * |_|  |_|\__,_|\__\___|_|  |_|\__,_|_| .__/ \___/|_|_| |_|\__|___/
   *                                     |_|
   * */

  std::unordered_map< std::string, MaterialPointFactory::materialPointFactoryFunction >& MaterialPointFactory::
    materialPointFactoryFunctionByName()
  {
    static std::unordered_map< std::string, materialPointFactoryFunction > map;
    return map;
  }

  bool MaterialPointFactory::registerMaterialPoint( const std::string&           materialPointName,
                                                          materialPointFactoryFunction factoryFunction )
  {

    const auto materialPointNameUpperCase = makeStringUpperCase( materialPointName );

    assert( materialPointFactoryFunctionByName().find( materialPointNameUpperCase ) ==
            materialPointFactoryFunctionByName().end() );

    materialPointFactoryFunctionByName()[materialPointNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Discretizations::MaterialPoints::MaterialPoint* MaterialPointFactory::createMaterialPoint( const std::string& materialPointName,
                                                                        int                materialPointNumber,
                                                                        const double*      vertexCoordinates,
                                                                        int                sizeVertexCoordinates,
                                                                        double             volume
                                                                        /* const Marmot::Discretizations::MaterialSection& material */
  )
  {
    const auto materialPointNameUpperCase = makeStringUpperCase( materialPointName );

    try {
      return materialPointFactoryFunctionByName().at(
        materialPointNameUpperCase )( materialPointNumber, vertexCoordinates, sizeVertexCoordinates, volume
                                      /* material */
      );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << "Invalid materialPoint " << materialPointName << " requested!" );
    }
  }

  /*
   *
   *  ____ _               _           _    ____     _ _
   * / ___| | __ _ ___ ___(_) ___ __ _| |  / ___|___| | |___
   *| |   | |/ _` / __/ __| |/ __/ _` | | | |   / _ \ | / __|
   *| |___| | (_| \__ \__ \ | (_| (_| | | | |__|  __/ | \__ \
   * \____|_|\__,_|___/___/_|\___\__,_|_|  \____\___|_|_|___/
   *
   */

  std::unordered_map< std::string, CellFactory::cellFactoryFunction >& CellFactory::
    cellFactoryFunctionByName()
  {
    static std::unordered_map< std::string, cellFactoryFunction > map;
    return map;
  }

  bool CellFactory::registerCell( const std::string& cellName, cellFactoryFunction factoryFunction )
  {

    const auto cellNameUpperCase = makeStringUpperCase( cellName );

    assert( cellFactoryFunctionByName().find( cellNameUpperCase ) == cellFactoryFunctionByName().end() );

    cellFactoryFunctionByName()[cellNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Discretizations::Cells::Cell* CellFactory::createCell( const std::string& cellName,
                                             int                cellNumber,
                                             const double*      nodeCoordinates,
                                             int                sizeNodeCoordinates )
  {
    const auto cellNameUpperCase = makeStringUpperCase( cellName );

    try {
      return cellFactoryFunctionByName().at( cellNameUpperCase )( cellNumber, nodeCoordinates, sizeNodeCoordinates );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << "Invalid cell " << cellName << " requested!" );
    }
  }

  /*
   * ____       ____        _ _               ____     _ _
   *| __ )     / ___| _ __ | (_)_ __   ___   / ___|___| | |___
   *|  _ \ ____\___ \| '_ \| | | '_ \ / _ \ | |   / _ \ | / __|
   *| |_) |_____|__) | |_) | | | | | |  __/ | |__|  __/ | \__ \
   *|____/     |____/| .__/|_|_|_| |_|\___|  \____\___|_|_|___/
   *                 |_|
   */

  std::unordered_map< std::string, CellFactory::bSplineCellFactoryFunction >& CellFactory::
    bSplineCellFactoryFunctionByName()
  {
    static std::unordered_map< std::string, bSplineCellFactoryFunction > map;
    return map;
  }

  bool CellFactory::registerBSplineCell( const std::string& cellName, bSplineCellFactoryFunction factoryFunction )
  {

    const auto cellNameUpperCase = makeStringUpperCase( cellName );

    assert( bSplineCellFactoryFunctionByName().find( cellNameUpperCase ) == bSplineCellFactoryFunctionByName().end() );

    bSplineCellFactoryFunctionByName()[cellNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Discretizations::Cells::Cell* CellFactory::createBSplineCell( const std::string& cellName,
                                                    int                cellNumber,
                                                    const double*      nodeCoordinates,
                                                    int                sizeNodeCoordinates,
                                                    const double*      knotVectors,
                                                    int                sizeKnotVectors )
  {
    const auto cellNameUpperCase = makeStringUpperCase( cellName );

    try {
      return bSplineCellFactoryFunctionByName().at(
        cellNameUpperCase )( cellNumber, nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << "Invalid cell " << cellName << " requested!" );
    }
  }

  /*
   *    ____     _ _ _____ _                           _
   *   / ___|___| | | ____| | ___ _ __ ___   ___ _ __ | |_ ___
   *  | |   / _ \ | |  _| | |/ _ \ '_ ` _ \ / _ \ '_ \| __/ __|
   *  | |__|  __/ | | |___| |  __/ | | | | |  __/ | | | |_\__ \
   *   \____\___|_|_|_____|_|\___|_| |_| |_|\___|_| |_|\__|___/
   */

  std::unordered_map< std::string, CellElementFactory::cellElementFactoryFunction >& CellElementFactory::
    cellElementFactoryFunctionByName()
  {
    static std::unordered_map< std::string, cellElementFactoryFunction > map;
    return map;
  }

  bool CellElementFactory::registerCellElement( const std::string&         cellElementName,
                                                      cellElementFactoryFunction factoryFunction )
  {

    const auto cellElementNameUpperCase = makeStringUpperCase( cellElementName );

    assert( cellElementFactoryFunctionByName().find( cellElementNameUpperCase ) ==
            cellElementFactoryFunctionByName().end() );

    cellElementFactoryFunctionByName()[cellElementNameUpperCase] = factoryFunction;

    return true;
  }

  Marmot::Discretizations::Cells::CellElement* CellElementFactory::createCellElement( const std::string& cellElementName,
                                                                  int                cellElementNumber,
                                                                  const double*      nodeCoordinates,
                                                                  int                sizeNodeCoordinates,
                                                                  const std::string& quadratureRule,
                                                                  int                quadratureOrder )
  {
    const auto cellNameUpperCase = makeStringUpperCase( cellElementName );

    try {
      return cellElementFactoryFunctionByName().at(
        cellNameUpperCase )( cellElementNumber, nodeCoordinates, sizeNodeCoordinates, quadratureRule, quadratureOrder );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << "Invalid cellElement " << cellElementName << " requested!" );
    }
  }

} // namespace Marmot::Registry
