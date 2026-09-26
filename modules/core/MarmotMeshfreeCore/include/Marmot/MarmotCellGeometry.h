/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
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

template < int nDim, int nNodes >
class [[deprecated( "Currently unused" )]] MarmotCellGeometry {

public:
  typedef Eigen::Matrix< double, nDim, 1 >          XiSized;
  typedef Eigen::Matrix< double, nDim * nNodes, 1 > CoordinateVector;
  typedef Eigen::Matrix< double, 1, nNodes >        NSized;
  typedef Eigen::Matrix< double, nDim, nNodes >     dNdXSized;

  virtual bool isCoordinateInCell( const double* coordinates ) const = 0;

  virtual XiSized findReferenceCoordinate( const XiSized& coord ) const = 0;

  virtual NSized N( const XiSized& xi ) const = 0;

  virtual dNdXSized dNdX( const XiSized& xi ) const = 0;

  virtual double detJ( const XiSized& xi ) const = 0;
};

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
