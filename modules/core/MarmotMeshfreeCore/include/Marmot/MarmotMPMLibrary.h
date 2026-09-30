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

  /**
   * @class MarmotLibrary::MarmotMaterialPointFactory
   * @brief Registry of MarmotMaterialPoint types, created by name.
   *
   * @details A material point module registers each of its types once, at static initialization time, with a
   * factory function (typically a captureless lambda) in its @c *Registration.cpp:
   * @code
   * const static bool isRegistered = MarmotLibrary::MarmotMaterialPointFactory::registerMaterialPoint(
   *   "Displacement/PlaneStrain",
   *   []( int n, const double* coords, int sizeCoords, double volume ) -> MarmotMaterialPoint* {
   *     return new DisplacementMaterialPoint2D( n, coords, sizeCoords, volume );
   *   } );
   * @endcode
   * Names are case-insensitive (stored in upper case). The class is not instantiable; all members are static.
   */
  class MarmotMaterialPointFactory {
  public:
    /**
     * @brief Signature of a material point factory function.
     * @details Arguments: material point number, vertex coordinates (vertex by vertex), size of the vertex
     * coordinate array, and the volume of the material point.
     */
    using materialPointFactoryFunction = MarmotMaterialPoint* (*)( int           materialPointNumber,
                                                                   const double* vertexCoordinates,
                                                                   int           sizeVertexCoordinates,
                                                                   double        volume
                                                                   /* const MarmotMaterialSection& material */
    );
    MarmotMaterialPointFactory()       = delete;

    /**
     * @brief Creates a material point of a registered type.
     * @param[in] materialPointName     Registered name (case-insensitive), e.g. @c "Displacement/PlaneStrain".
     * @param[in] materialNumber        Number of the new material point.
     * @param[in] vertexCoordinates     Vertex coordinates (vertex by vertex).
     * @param[in] sizeVertexCoordinates Size of @p vertexCoordinates.
     * @param[in] volume                Volume of the material point.
     * @return Newly allocated material point; the caller takes ownership.
     * @throws std::invalid_argument if no material point is registered under @p materialPointName.
     */
    static MarmotMaterialPoint* createMaterialPoint( const std::string& materialPointName,
                                                     int                materialNumber,
                                                     const double*      vertexCoordinates,
                                                     int                sizeVertexCoordinates,
                                                     double             volume
                                                     /* const MarmotMaterialSection& material */
    );

    /**
     * @brief Registers a material point type.
     * @param[in] materialName    Name (case-insensitive); registering a name twice fails an @c assert.
     * @param[in] factoryFunction Function creating an instance.
     * @return @c true (used to initialize a static flag).
     */
    static bool registerMaterialPoint( const std::string& materialName, materialPointFactoryFunction factoryFunction );

  private:
    /**
     * @brief Checks whether a name is registered (declared, but not implemented).
     * @param[in] materialPointName Name of the material point type.
     * @return @c true if registered.
     */
    bool checkIfMaterialPointIsRegistered( const std::string& materialPointName );

    /// Factory functions by upper-case name (a function-local static: registrations run during static initialization,
    /// possibly before a static data member of this translation unit would be constructed).
    static std::unordered_map< std::string, materialPointFactoryFunction >& materialPointFactoryFunctionByName();
  };

  /**
   * @class MarmotLibrary::MarmotCellFactory
   * @brief Registry of MarmotCell types, created by name.
   *
   * @details Keeps two separate registries: Lagrangian cells, created from their node coordinates
   * (registerCell() / createCell(), e.g. @c "Displacement/Quad4"), and B-spline cells, which additionally need
   * the knot vectors of the cell (registerBSplineCell() / createBSplineCell(), e.g. @c "Displacement/BSpline/2").
   * Registration works as for MarmotMaterialPointFactory. Names are case-insensitive (stored in upper case).
   */
  class MarmotCellFactory {
  public:
    /**
     * @brief Signature of a (Lagrangian) cell factory function.
     * @details Arguments: cell number, node coordinates (node by node), size of the coordinate array.
     */
    using cellFactoryFunction = MarmotCell* (*)( int           cellNumber,
                                                 const double* vertexCoordinates,
                                                 int           sizeVertexCoordinates );

    /**
     * @brief Signature of a B-spline cell factory function.
     * @details Arguments: cell number, control point coordinates (point by point), size of the coordinate array,
     * knot vectors (direction by direction), size of the knot vector array.
     */
    using bSplineCellFactoryFunction = MarmotCell* (*)( int           cellNumber,
                                                        const double* vertexCoordinates,
                                                        int           sizeVertexCoordinates,
                                                        const double* knotVectors,
                                                        int           sizeKnotVectors );

    MarmotCellFactory() = delete;

    /**
     * @brief Creates a Lagrangian cell of a registered type.
     * @param[in] cellName            Registered name (case-insensitive), e.g. @c "Displacement/Quad4".
     * @param[in] cellNumber          Number of the new cell.
     * @param[in] nodeCoordinates     Node coordinates (node by node); the cells map this array, so it must
     *                                outlive the cell.
     * @param[in] sizeNodeCoordinates Size of @p nodeCoordinates.
     * @return Newly allocated cell; the caller takes ownership.
     * @throws std::invalid_argument if no cell is registered under @p cellName.
     */
    static MarmotCell* createCell( const std::string& cellName,
                                   int                cellNumber,
                                   const double*      nodeCoordinates,
                                   int                sizeNodeCoordinates );

    /**
     * @brief Creates a B-spline cell of a registered type.
     * @param[in] cellName            Registered name (case-insensitive), e.g. @c "Displacement/BSpline/2".
     * @param[in] cellNumber          Number of the new cell.
     * @param[in] nodeCoordinates     Control point coordinates (point by point); mapped by the cell, so it must
     *                                outlive the cell.
     * @param[in] sizeNodeCoordinates Size of @p nodeCoordinates.
     * @param[in] knotVectors         Knot vectors of the cell, @f$ 2p+2 @f$ knots per direction, direction by
     *                                direction.
     * @param[in] sizeKnotVectors     Size of @p knotVectors.
     * @return Newly allocated cell; the caller takes ownership.
     * @throws std::invalid_argument if no B-spline cell is registered under @p cellName.
     */
    static MarmotCell* createBSplineCell( const std::string& cellName,
                                          int                cellNumber,
                                          const double*      nodeCoordinates,
                                          int                sizeNodeCoordinates,
                                          const double*      knotVectors,
                                          int                sizeKnotVectors );

    /**
     * @brief Registers a Lagrangian cell type.
     * @param[in] cellName        Name (case-insensitive); registering a name twice fails an @c assert.
     * @param[in] factoryFunction Function creating an instance.
     * @return @c true (used to initialize a static flag).
     */
    static bool registerCell( const std::string& cellName, cellFactoryFunction factoryFunction );
    /**
     * @brief Registers a B-spline cell type.
     * @param[in] cellName        Name (case-insensitive); registering a B-spline name twice fails an @c assert, a
     *                            Lagrangian cell of the same name is independent.
     * @param[in] factoryFunction Function creating an instance.
     * @return @c true (used to initialize a static flag).
     */
    static bool registerBSplineCell( const std::string& cellName, bSplineCellFactoryFunction factoryFunction );

  private:
    /// Lagrangian cell factory functions by upper-case name (a function-local static: registrations run during static
    /// initialization, possibly before a static data member of this translation unit would be constructed).
    static std::unordered_map< std::string, cellFactoryFunction >& cellFactoryFunctionByName();
    /// B-spline cell factory functions by upper-case name (a function-local static, see above).
    static std::unordered_map< std::string, bSplineCellFactoryFunction >& bSplineCellFactoryFunctionByName();
  };

  /**
   * @class MarmotLibrary::MarmotCellElementFactory
   * @brief Registry of MarmotCellElement types, created by name.
   *
   * @details Registration works as for MarmotMaterialPointFactory; in addition to the node coordinates, a cell
   * element receives the quadrature rule that defines its requested material points. Names are case-insensitive
   * (stored in upper case). No cell element is registered in Marmot at present.
   */
  class MarmotCellElementFactory {
  public:
    /**
     * @brief Signature of a cell element factory function.
     * @details Arguments: cell element number, node coordinates (node by node), size of the coordinate array,
     * name of the quadrature rule, quadrature order.
     */
    using cellElementFactoryFunction = MarmotCellElement* (*)( int                cellElementNumber,
                                                               const double*      vertexCoordinates,
                                                               int                sizeVertexCoordinates,
                                                               const std::string& quadratureRule,
                                                               int                quadratureOrder );

    MarmotCellElementFactory() = delete;

    /**
     * @brief Creates a cell element of a registered type.
     * @param[in] cellElementName     Registered name (case-insensitive).
     * @param[in] cellElementNumber   Number of the new cell element.
     * @param[in] nodeCoordinates     Node coordinates (node by node).
     * @param[in] sizeNodeCoordinates Size of @p nodeCoordinates.
     * @param[in] quadratureRule      Name of the quadrature rule for the requested material points.
     * @param[in] quadratureOrder     Order of the quadrature rule.
     * @return Newly allocated cell element; the caller takes ownership.
     * @throws std::invalid_argument if no cell element is registered under @p cellElementName.
     */
    static MarmotCellElement* createCellElement( const std::string& cellElementName,
                                                 int                cellElementNumber,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 const std::string& quadratureRule,
                                                 int                quadratureOrder );

    /**
     * @brief Registers a cell element type.
     * @param[in] cellElementName Name (case-insensitive); registering a name twice fails an @c assert.
     * @param[in] factoryFunction Function creating an instance.
     * @return @c true (used to initialize a static flag).
     */
    static bool registerCellElement( const std::string& cellElementName, cellElementFactoryFunction factoryFunction );

  private:
    /// Cell element factory functions by upper-case name (a function-local static: registrations run during static
    /// initialization, possibly before a static data member of this translation unit would be constructed).
    static std::unordered_map< std::string, cellElementFactoryFunction >& cellElementFactoryFunctionByName();
  };

} // namespace MarmotLibrary
