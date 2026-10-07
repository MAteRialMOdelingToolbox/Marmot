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
#include "Marmot/MarmotBSplineGeometryElement.h"
#include "Marmot/MarmotCellGeometry.h"
#include <iostream>

namespace Marmot::Cells {

  /**
   * @class Marmot::Cells::MarmotBSplineCellGeometry
   * @brief B-spline cell geometry (one knot span), a geometry policy for the MPM cells.
   *
   * @details Wraps Marmot::FiniteElement::MarmotBSplineGeometryElement for use as the @c GeometryCell of a cell
   * (see GeometryCellPolicy). The cell is the knot span
   * @f$ [u_{p}, u_{p+1}] @f$ of each direction; the knot vectors are given in **physical coordinates**, so the
   * parametric coordinate of a point is the point itself: findReferenceCoordinate() is the identity and dNdX()
   * returns the parametric derivatives unchanged. This is exact for a B-spline grid whose geometry map is the
   * identity @f$ \boldsymbol{X}(\boldsymbol{\xi}) = \boldsymbol{\xi} @f$ (control points at the Greville abscissae,
   * e.g. a uniform, axis-aligned background grid); the control point coordinates enter only detJ().
   *
   * @tparam nDim  Spatial dimension (2 or 3).
   * @tparam order Polynomial degree @f$ p @f$ (1, 2 or 3 are registered).
   */
  template < int nDim, int order >
  class MarmotBSplineCellGeometry : public Marmot::FiniteElement::MarmotBSplineGeometryElement< nDim, order >

  {

    /// The underlying B-spline geometry element.
    using ParentBSplineGeometryElement = Marmot::FiniteElement::MarmotBSplineGeometryElement< nDim, order >;
    using JacobianSized                = ParentBSplineGeometryElement::JacobianSized; ///< Jacobian matrix type.

    Eigen::Matrix< double, nDim, 1 > _boundingBoxMin; ///< Lower corner of the knot span, @f$ u_p @f$ per direction.
    Eigen::Matrix< double, nDim, 1 > _boundingBoxMax; ///< Upper corner of the knot span, @f$ u_{p+1} @f$ per direction.

    bool _boundingBoxMatchesGeometryExactly;          ///< Always @c true; currently not used.

  public:
    using NSized           = ParentBSplineGeometryElement::NSized;           ///< Row vector of shape function values.
    using dNdXSized        = ParentBSplineGeometryElement::dNdXSized;        ///< Shape function gradients.
    using CoordinateVector = ParentBSplineGeometryElement::CoordinateVector; ///< Flat control point coordinates.
    using XiSized          = ParentBSplineGeometryElement::XiSized;          ///< Parametric (= physical) point.

    /**
     * @brief Constructs the geometry and its bounding box (the knot span).
     * @param[in] nodeCoordinates     Control point coordinates, point by point; mapped, not copied.
     * @param[in] sizeNodeCoordinates Size of @p nodeCoordinates.
     * @param[in] knotVectors         Knot vectors (@f$ 2p+2 @f$ knots per direction, direction by direction), in
     *                                physical coordinates.
     * @param[in] sizeKnotVectors     Size of @p knotVectors.
     */
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

    /**
     * @brief Knot span test @f$ u_{p,i} \le x_i < u_{p+1,i} @f$ (half-open).
     * @param[in] coordinates Point coordinates (nDim values).
     * @return @c true if the point is inside the knot span.
     */
    bool isCoordinateInCell( const double* coordinates ) const;

    /**
     * @brief Bounding box of the cell, i.e. the knot span.
     * @param[out] boundingBoxMin Lower corner (nDim values).
     * @param[out] boundingBoxMax Upper corner (nDim values).
     */
    void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const;

    /**
     * @brief Inverse geometry map; the identity, since the knots are physical coordinates.
     * @param[in] coord Physical coordinates.
     * @return @p coord.
     */
    XiSized findReferenceCoordinate( const XiSized& coord ) const;

    /**
     * @brief Tensor-product B-spline shape functions.
     * @param[in] xi Parametric (= physical) coordinates.
     * @return @f$ N_a(\boldsymbol{\xi}) @f$.
     */
    NSized N( const XiSized& xi ) const { return ParentBSplineGeometryElement::N( xi ); }

    /**
     * @brief Shape function gradients; the parametric derivatives, no Jacobian transformation (identity map).
     * @param[in] xi Parametric (= physical) coordinates.
     * @return @f$ \partial N_a/\partial \xi_i @f$ (row @f$ i @f$, column @f$ a @f$).
     */
    dNdXSized dNdX( const XiSized& xi ) const
    {
      const auto dN_dXi = ParentBSplineGeometryElement::dNdXi( xi );
      return dN_dXi;
    }

    /**
     * @brief Determinant of the Jacobian of the geometry map defined by the control points.
     * @param[in] xi Parametric coordinates.
     * @return @f$ \det \boldsymbol{J}(\boldsymbol{\xi}) @f$.
     */
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
