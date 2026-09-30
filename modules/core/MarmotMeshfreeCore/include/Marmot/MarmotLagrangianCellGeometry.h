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

/**
 * @class MarmotLagrangianCellGeometry
 * @brief Isoparametric Lagrangian cell geometry (Quad4, Hexa8), a geometry policy for the MPM cells.
 *
 * @details Uses the shape functions @f$ N_A(\boldsymbol{\xi}) @f$, @f$ \boldsymbol{\xi}\in[-1,1]^{n_\mathrm{dim}} @f$,
 * of MarmotGeometryElement and the isoparametric map
 * @f$ \boldsymbol{X}(\boldsymbol{\xi}) = \sum_A N_A(\boldsymbol{\xi})\,\boldsymbol{X}_A @f$ with the nodal coordinates
 * @f$ \boldsymbol{X}_A @f$, and provides the inverse map, the physical gradients and the point location test
 * required by GeometryCellPolicy.
 *
 * @warning The inverse map findReferenceCoordinate() and the point location test isCoordinateInCell() are based on
 * the axis-aligned bounding box of the nodes, and are therefore only valid for **axis-aligned rectangular
 * (box-shaped) cells**. For a distorted cell, findReferenceCoordinate() throws (the Newton update of the inverse map is
 * not implemented yet).
 *
 * @tparam nDim   Spatial dimension (2 or 3).
 * @tparam nNodes Number of nodes (4 in 2D, 8 in 3D).
 */
template < int nDim, int nNodes >
class MarmotLagrangianCellGeometry : public MarmotGeometryElement< nDim, nNodes >

{

  using ParentLagrangianGeometryElement = MarmotGeometryElement< nDim, nNodes >; ///< Lagrangian geometry element.
  using JacobianSized = typename ParentLagrangianGeometryElement::JacobianSized; ///< Jacobian matrix type.
  Eigen::Matrix< double, nDim, 1 > _boundingBoxMin;                              ///< Lower corner of the bounding box.
  Eigen::Matrix< double, nDim, 1 > _boundingBoxMax;                              ///< Upper corner of the bounding box.

  bool _boundingBoxMatchesGeometryExactly; ///< Always @c true (assumes a box-shaped cell); currently not used.

public:
  using NSized    = typename ParentLagrangianGeometryElement::NSized;     ///< Row vector of shape function values.
  using dNdXSized = typename ParentLagrangianGeometryElement::dNdXiSized; ///< Shape function gradients (nDim x nNodes).
  using XiSized   = typename ParentLagrangianGeometryElement::XiSized;    ///< Parametric (or physical) point.

  /**
   * @brief Constructs the geometry and its bounding box.
   * @param[in] nodeCoordinates Nodal coordinates, nNodes @f$ \times @f$ nDim values, node by node (Abaqus node
   *                            order). The array is mapped, not copied, and must outlive the geometry.
   */
  MarmotLagrangianCellGeometry( const double* nodeCoordinates )
  {
    ParentLagrangianGeometryElement::assignNodeCoordinates( nodeCoordinates );
    auto nodeCoords = ParentLagrangianGeometryElement::coordinates.reshaped( nDim, nNodes );

    _boundingBoxMin = nodeCoords.rowwise().minCoeff();
    _boundingBoxMax = nodeCoords.rowwise().maxCoeff();

    _boundingBoxMatchesGeometryExactly = true;
  }

  /**
   * @brief Bounding box test @f$ X_{\min,i} \le x_i < X_{\max,i} @f$ (half-open, so that a point on a shared face
   * belongs to exactly one of two neighbouring cells). Exact only for axis-aligned box-shaped cells.
   * @param[in] coordinates Point coordinates (nDim values).
   * @return @c true if the point is inside the bounding box.
   */
  bool isCoordinateInCell( const double* coordinates ) const;

  /**
   * @brief Axis-aligned bounding box of the nodes.
   * @param[out] boundingBoxMin Lower corner (nDim values).
   * @param[out] boundingBoxMax Upper corner (nDim values).
   */
  void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const;

  /**
   * @brief Inverse isoparametric map @f$ \boldsymbol{\xi}(\boldsymbol{X}) @f$.
   * @details Uses the affine map of the bounding box,
   * @f$ \xi_i = 2\,(X_i - \tfrac12(X_{\max,i}+X_{\min,i}))/(X_{\max,i}-X_{\min,i}) @f$, and checks the residual
   * @f$ \|\boldsymbol{X} - \sum_A N_A(\boldsymbol{\xi})\boldsymbol{X}_A\| / \|\boldsymbol{X}\| < 10^{-12} @f$. There is
   * no Newton update yet, so the result is correct only for axis-aligned box-shaped cells.
   * @param[in] coord Physical coordinates.
   * @return Parametric coordinates @f$ \boldsymbol{\xi} @f$.
   * @throws std::runtime_error if the residual check fails (distorted cell).
   */
  XiSized findReferenceCoordinate( const XiSized& coord ) const;

  /**
   * @brief Shape functions.
   * @param[in] xi Parametric coordinates.
   * @return @f$ N_A(\boldsymbol{\xi}) @f$.
   */
  NSized N( const XiSized& xi ) const { return ParentLagrangianGeometryElement::N( xi ); }

  /**
   * @brief Shape function gradients w.r.t. the physical coordinates,
   * @f$ \partial N_A/\partial X_i = (\partial N_A/\partial \xi_j)\,J^{-1}_{ji} @f$,
   * @f$ J_{ij} = \partial X_i/\partial \xi_j @f$.
   * @param[in] xi Parametric coordinates.
   * @return @f$ \partial N_A/\partial X_i @f$ (row @f$ i @f$, column @f$ A @f$).
   */
  dNdXSized dNdX( const XiSized& xi ) const
  {
    const auto          dN_dXi = ParentLagrangianGeometryElement::dNdXi( xi );
    const JacobianSized J      = ParentLagrangianGeometryElement::Jacobian( dN_dXi );
    const JacobianSized invJ   = J.inverse();

    return ParentLagrangianGeometryElement::dNdX( dN_dXi, invJ );
  }

  /**
   * @brief Determinant of the Jacobian of the isoparametric map.
   * @param[in] xi Parametric coordinates.
   * @return @f$ \det \boldsymbol{J}(\boldsymbol{\xi}) @f$.
   */
  double detJ( const XiSized& xi ) const
  {
    const auto dN_dXi = ParentLagrangianGeometryElement::dNdXi( xi );
    return ParentLagrangianGeometryElement::Jacobian( dN_dXi ).determinant();
  };

  /**
   * @brief Placeholder (unused).
   * @param[in] xi Parametric coordinates (ignored).
   * @return Always @c true.
   */
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
