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
#include "Marmot/MarmotBSplineGeometryElement.h"
#include "Marmot/MarmotCellGeometry.h"
#include <iostream>

namespace Marmot::Cells {

  template < int nDim, int order >
  class MarmotBSplineCellGeometry : public Marmot::FiniteElement::MarmotBSplineGeometryElement< nDim, order >

  {

    using ParentBSplineGeometryElement = Marmot::FiniteElement::MarmotBSplineGeometryElement< nDim, order >;
    using JacobianSized                = ParentBSplineGeometryElement::JacobianSized;

    Eigen::Matrix< double, nDim, 1 > _boundingBoxMin;
    Eigen::Matrix< double, nDim, 1 > _boundingBoxMax;

    bool _boundingBoxMatchesGeometryExactly;

  public:
    using NSized           = ParentBSplineGeometryElement::NSized;
    using dNdXSized        = ParentBSplineGeometryElement::dNdXSized;
    using CoordinateVector = ParentBSplineGeometryElement::CoordinateVector;
    using XiSized          = ParentBSplineGeometryElement::XiSized;

    MarmotBSplineCellGeometry( const double* nodeCoordinates,
                               int           sizeNodeCoordinates,
                               const double* knotVectors,
                               int           sizeKnotVectors )
      : ParentBSplineGeometryElement( nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors )
    {
      auto nodeCoords = ParentBSplineGeometryElement::_mapCoordinates.reshaped( nDim, this->nNodes );

      /* _boundingBoxMin = nodeCoords.rowwise().minCoeff(); */
      /* _boundingBoxMax = nodeCoords.rowwise().maxCoeff(); */
      _boundingBoxMin = this->_knotVectors.row( order );
      _boundingBoxMax = this->_knotVectors.row( this->nKnotsPerDir - order - 1 );

      _boundingBoxMatchesGeometryExactly = true;
    }

    bool isCoordinateInCell( const double* coordinates ) const;

    void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const;

    XiSized findReferenceCoordinate( const XiSized& coord ) const;

    NSized N( const XiSized& xi ) const { return ParentBSplineGeometryElement::N( xi ); }

    dNdXSized dNdX( const XiSized& xi ) const
    {
      const auto dN_dXi = ParentBSplineGeometryElement::dNdXi( xi );
      return dN_dXi;
    }

    double detJ( const XiSized& xi ) const
    {
      const auto dN_dXi = ParentBSplineGeometryElement::dNdXi( xi );
      return ParentBSplineGeometryElement::Jacobian( dN_dXi ).determinant();
    };
  };

  template < int nDim, int order >
  bool MarmotBSplineCellGeometry< nDim, order >::isCoordinateInCell( const double* coordinates ) const
  {

    for ( auto i = 0; i < nDim; i++ )
      if ( coordinates[i] < _boundingBoxMin( i ) || coordinates[i] >= _boundingBoxMax( i ) )
        return false;

    return true;
  }

  template < int nDim, int nNodes >
  void MarmotBSplineCellGeometry< nDim, nNodes >::getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const
  {
    ( Eigen::Map< XiSized >( boundingBoxMin ) ) = _boundingBoxMin;
    ( Eigen::Map< XiSized >( boundingBoxMax ) ) = _boundingBoxMax;
  }

  template < int nDim, int order >
  MarmotBSplineCellGeometry< nDim, order >::XiSized MarmotBSplineCellGeometry< nDim, order >::findReferenceCoordinate(
    const XiSized& coord ) const
  {
    // initial guess:
    /* XiSized xi = 2 * ( coord - ( _boundingBoxMax + _boundingBoxMin ) / 2 ) */
    /*                    .cwiseProduct( ( _boundingBoxMax - _boundingBoxMin ).cwiseInverse() ); */

    return coord;

    /*     std::cout << "coord " << coord.transpose() << std::endl; */
    /*     std::cout << "initial guess " << xi.transpose() << std::endl; */
    /*     std::cout << "bb " << _boundingBoxMin.transpose() << std::endl; */
    /*     std::cout << "bb " << _boundingBoxMax.transpose() << std::endl; */

    /* XiSized r = coord - ParentBSplineGeometryElement::_mapCoordinates.reshaped(nDim, this->nNodes) * N( xi
     * ).transpose(); */

    /* int nCounter = 0; */
    /* while ( r.norm() / coord.norm() >= 1e-12 ) { */

    /*   // TODO */
    /*   /1* xi += *1/ */

    /*   r = coord - ParentBSplineGeometryElement::_mapCoordinates.reshaped(nDim, this->nNodes) * N( xi ).transpose();
     */
    /*   nCounter++; */
    /*   if ( nCounter >= 5 ) { */
    /*     throw std::runtime_error( MakeString() */
    /*                               << __PRETTY_FUNCTION__ << ": failed to determine inverse map for coordinate " */
    /*                               << coord.transpose() ); */
    /*   } */
    /* } */

    /* return xi; */
  }

} // namespace Marmot::Cells
