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
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include <concepts>
#include <stdexcept>

/**
 * @class MarmotCellGeometry
 * @brief Deprecated dynamic-polymorphism interface of a cell geometry; superseded by GeometryCellPolicy.
 *
 * @details Unused: the cells take their geometry as a template parameter constrained by the concept
 * GeometryCellPolicy (static polymorphism), with MarmotLagrangianCellGeometry and
 * Marmot::Cells::MarmotBSplineCellGeometry as implementations. The members describe the same contract.
 *
 * @tparam nDim   Spatial dimension.
 * @tparam nNodes Number of cell nodes.
 */
template < int nDim, int nNodes >
class [[deprecated( "Currently unused" )]] MarmotCellGeometry {

public:
  typedef Eigen::Matrix< double, nDim, 1 >          XiSized;          ///< Parametric (or physical) point.
  typedef Eigen::Matrix< double, nDim * nNodes, 1 > CoordinateVector; ///< Flat nodal coordinates.
  typedef Eigen::Matrix< double, 1, nNodes >        NSized;           ///< Row vector of shape function values.
  typedef Eigen::Matrix< double, nDim, nNodes >     dNdXSized;        ///< Shape function gradients.

  /**
   * @brief Checks whether a point lies inside the cell.
   * @param[in] coordinates Point coordinates (nDim values).
   * @return @c true if the point is inside the cell.
   */
  virtual bool isCoordinateInCell( const double* coordinates ) const = 0;

  /**
   * @brief Inverse geometry map: parametric coordinates of a physical point.
   * @param[in] coord Physical coordinates.
   * @return Parametric coordinates @f$ \boldsymbol{\xi} @f$.
   */
  virtual XiSized findReferenceCoordinate( const XiSized& coord ) const = 0;

  /**
   * @brief Shape functions.
   * @param[in] xi Parametric coordinates.
   * @return @f$ N_A(\boldsymbol{\xi}) @f$.
   */
  virtual NSized N( const XiSized& xi ) const = 0;

  /**
   * @brief Shape function gradients w.r.t. the physical coordinates.
   * @param[in] xi Parametric coordinates.
   * @return @f$ \partial N_A / \partial X_i @f$ (row @f$ i @f$, column @f$ A @f$).
   */
  virtual dNdXSized dNdX( const XiSized& xi ) const = 0;

  /**
   * @brief Determinant of the Jacobian of the geometry map.
   * @param[in] xi Parametric coordinates.
   * @return @f$ \det(\partial \boldsymbol{X}/\partial\boldsymbol{\xi}) @f$.
   */
  virtual double detJ( const XiSized& xi ) const = 0;
};

/**
 * @brief Requirements on the geometry policy of a cell (e.g. Marmot::Cells::DisplacementCell).
 *
 * @details A geometry policy must provide the types @c XiSized (@f$ n_\mathrm{dim} @f$ vector), @c NSized
 * (@f$ 1\times n_\mathrm{nodes} @f$) and @c dNdXSized (@f$ n_\mathrm{dim}\times n_\mathrm{nodes} @f$) and the
 * methods @c findReferenceCoordinate, @c N, @c dNdX and @c detJ with the signatures of MarmotCellGeometry.
 * The cells additionally call @c isCoordinateInCell, @c getBoundingBox and @c getElementShape, which the
 * concept does not check. Implementations: MarmotLagrangianCellGeometry, Marmot::Cells::MarmotBSplineCellGeometry.
 *
 * @tparam GeometryCellImpl Candidate geometry class.
 * @tparam nDim             Spatial dimension.
 * @tparam nNodes           Number of cell nodes.
 */
template < class GeometryCellImpl, int nDim, int nNodes >
concept GeometryCellPolicy = requires( GeometryCellImpl geom ) {
  typename GeometryCellImpl::XiSized;
  requires std::same_as< typename GeometryCellImpl::XiSized, typename Eigen::Matrix< double, nDim, 1 > >;

  typename GeometryCellImpl::NSized;
  requires std::same_as< typename GeometryCellImpl::NSized, typename Eigen::Matrix< double, 1, nNodes > >;

  typename GeometryCellImpl::dNdXSized;
  requires std::same_as< typename GeometryCellImpl::dNdXSized, typename Eigen::Matrix< double, nDim, nNodes > >;

  /* { */
  /*   geom.test( typename GeometryCellImpl::XiSized() ) */
  /* } -> std::same_as< bool >; */

  /* { */
  /*   geom.isCoordinateInCell( typename GeometryCellImpl::XiSized() ) */
  /* } -> std::same_as< bool >; */

  {
    geom.findReferenceCoordinate( typename GeometryCellImpl::XiSized() )
  } -> std::same_as< typename GeometryCellImpl::XiSized >;

  {
    geom.N( typename GeometryCellImpl::XiSized() )
  } -> std::same_as< typename GeometryCellImpl::NSized >;

  {
    geom.dNdX( typename GeometryCellImpl::XiSized() )
  } -> std::same_as< typename GeometryCellImpl::dNdXSized >;

  {
    geom.detJ( typename GeometryCellImpl::XiSized() )
  } -> std::same_as< double >;
};
