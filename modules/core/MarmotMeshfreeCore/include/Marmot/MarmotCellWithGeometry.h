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
#pragma once
#include "Marmot/MarmotCell.h"
#include <Eigen/Core>

namespace Marmot::Discretizations::Cells {

  /**
   * @class Marmot::Discretizations::Cells::CellWithGeometry
   * @brief Binds a geometry policy to the Cell interface and implements the purely geometric members of the
   * interface, independently of the physics.
   *
   * @details A physics cell derives from CellWithGeometry< CellBase, GeometryCell > instead of from CellBase and
   * GeometryCell, and inherits getInterpolationVector() from the geometry (its NAtCoordinates()), so
   * that every cell evaluates its basis the same way, inside the cell and, extrapolated, outside of it (the basis
   * extrapolation of the aggregated MPM, Coombs et al., CMAME 442 (2025) 118012). A cell may still override it.
   *
   * @tparam CellBase     Base class, Marmot::Discretizations::Cells::Cell.
   * @tparam GeometryCell Geometry policy (see GeometryCellPolicy), Lagrangian or B-spline.
   */
  template < class CellBase, class GeometryCell >
  class CellWithGeometry : public CellBase, public GeometryCell {

  public:
    /**
     * @brief Constructs the cell base and copies the geometry.
     * @param[in] geometry Geometry of the cell.
     */
    explicit CellWithGeometry( const GeometryCell& geometry ) : GeometryCell( geometry ){};

    /**
     * @brief Shape functions of the cell at an arbitrary point, also outside of the cell: there, the polynomial
     * extension of the cell's basis (no clamping), see NAtCoordinates() of the geometry policies.
     * @param[out] vec Shape functions (getNNodes() values).
     * @param[in] coordinates Coordinates of the point (nDim values).
     */
    void getInterpolationVector( double* vec, const double* coordinates ) const override
    {
      using XiSized = typename GeometryCell::XiSized;
      using NSized  = typename GeometryCell::NSized;
      Eigen::Map< NSized > interpolationVector( vec );
      interpolationVector = GeometryCell::NAtCoordinates( Eigen::Map< const XiSized >( coordinates ) );
    }
  };

} // namespace Marmot::Discretizations::Cells
