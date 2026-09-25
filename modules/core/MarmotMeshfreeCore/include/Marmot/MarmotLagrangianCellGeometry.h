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
#include "Marmot/MarmotCellGeometry.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include <stdexcept>

template < int nDim, int nNodes >
class MarmotLagrangianCellGeometry : public MarmotGeometryElement< nDim, nNodes >

{

  using ParentLagrangianGeometryElement = MarmotGeometryElement< nDim, nNodes >;
  using JacobianSized                   = typename ParentLagrangianGeometryElement::JacobianSized;
  Eigen::Matrix< double, nDim, 1 > _boundingBoxMin;
  Eigen::Matrix< double, nDim, 1 > _boundingBoxMax;

  bool _boundingBoxMatchesGeometryExactly;

public:
  using NSized    = typename ParentLagrangianGeometryElement::NSized;
  using dNdXSized = typename ParentLagrangianGeometryElement::dNdXiSized;
  using XiSized   = typename ParentLagrangianGeometryElement::XiSized;

  MarmotLagrangianCellGeometry( const double* nodeCoordinates )
  {
    ParentLagrangianGeometryElement::assignNodeCoordinates( nodeCoordinates );
    auto nodeCoords = ParentLagrangianGeometryElement::coordinates.reshaped( nDim, nNodes );

    _boundingBoxMin = nodeCoords.rowwise().minCoeff();
    _boundingBoxMax = nodeCoords.rowwise().maxCoeff();

    _boundingBoxMatchesGeometryExactly = true;
  }

  bool isCoordinateInCell( const double* coordinates ) const;

  void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const;

  XiSized findReferenceCoordinate( const XiSized& coord ) const;

  NSized N( const XiSized& xi ) const { return ParentLagrangianGeometryElement::N( xi ); }

  dNdXSized dNdX( const XiSized& xi ) const
  {
    const auto          dN_dXi = ParentLagrangianGeometryElement::dNdXi( xi );
    const JacobianSized J      = ParentLagrangianGeometryElement::Jacobian( dN_dXi );
    const JacobianSized invJ   = J.inverse();

    return ParentLagrangianGeometryElement::dNdX( dN_dXi, invJ );
  }

  double detJ( const XiSized& xi ) const
  {
    const auto dN_dXi = ParentLagrangianGeometryElement::dNdXi( xi );
    return ParentLagrangianGeometryElement::Jacobian( dN_dXi ).determinant();
  };

  bool test( const XiSized& xi ) { return true; };
};

template < int nDim, int nNodes >
bool MarmotLagrangianCellGeometry< nDim, nNodes >::isCoordinateInCell( const double* coordinates ) const
{

  for ( auto i = 0; i < nDim; i++ )
    if ( coordinates[i] < _boundingBoxMin( i ) || coordinates[i] >= _boundingBoxMax( i ) )
      return false;

  return true;
}

template < int nDim, int nNodes >
void MarmotLagrangianCellGeometry< nDim, nNodes >::getBoundingBox( double* boundingBoxMin,
                                                                   double* boundingBoxMax ) const
{
  ( Eigen::Map< XiSized >( boundingBoxMin ) ) = _boundingBoxMin;
  ( Eigen::Map< XiSized >( boundingBoxMax ) ) = _boundingBoxMax;
}

template < int nDim, int nNodes >
MarmotLagrangianCellGeometry< nDim, nNodes >::XiSized MarmotLagrangianCellGeometry< nDim, nNodes >::
  findReferenceCoordinate( const XiSized& coord ) const
{
  // initial guess:
  XiSized xi = 2 * ( coord - ( _boundingBoxMax + _boundingBoxMin ) / 2 )
                     .cwiseProduct( ( _boundingBoxMax - _boundingBoxMin ).cwiseInverse() );

  XiSized r = coord - ParentLagrangianGeometryElement::coordinates.reshaped( nDim, nNodes ) * N( xi ).transpose();

  int nCounter = 0;
  while ( r.norm() / coord.norm() >= 1e-12 ) {

    // TODO
    /* xi += */

    r = coord - ParentLagrangianGeometryElement::coordinates.reshaped( nDim, nNodes ) * N( xi ).transpose();
    nCounter++;
    if ( nCounter >= 5 ) {
      throw std::runtime_error( MakeString()
                                << __PRETTY_FUNCTION__ << ": failed to determine inverse map for coordinate "
                                << coord.transpose() );
    }
  }

  return xi;
}
