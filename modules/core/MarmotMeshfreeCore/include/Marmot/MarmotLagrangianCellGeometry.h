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
#include "Marmot/MarmotCellGeometry.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include <cmath>
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
 * For axis-aligned rectangular (box-shaped) cells, the inverse map findReferenceCoordinate() and the point location
 * test isCoordinateInCell() are exact and cheap (affine map of the bounding box); for distorted cells, the inverse map
 * is computed by Newton's method, and the point location test uses it.
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

  bool _boundingBoxMatchesGeometryExactly; ///< @c true for an axis-aligned box-shaped cell (every node is a box corner)

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

    // a box-shaped cell has all its nodes at corners of its bounding box
    const double tol                   = 1e-12 * ( _boundingBoxMax - _boundingBoxMin ).norm();
    _boundingBoxMatchesGeometryExactly = true;
    for ( int A = 0; A < nNodes; A++ )
      for ( int i = 0; i < nDim; i++ )
        if ( std::abs( nodeCoords( i, A ) - _boundingBoxMin( i ) ) > tol &&
             std::abs( nodeCoords( i, A ) - _boundingBoxMax( i ) ) > tol )
          _boundingBoxMatchesGeometryExactly = false;
  }

  /**
   * @brief Point location test: bounding box test @f$ X_{\min,i} \le x_i < X_{\max,i} @f$, and for a distorted
   * cell additionally @f$ -1 \le \xi_i < 1 @f$ for the parametric coordinates of findReferenceCoordinate().
   * @details Both tests are half-open, so that a point on a face shared by two cells belongs to exactly one of them;
   * consequently, a point exactly on the upper boundary of the whole mesh belongs to no cell.
   * @param[in] coordinates Point coordinates (nDim values).
   * @return @c true if the point is inside the cell.
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
   * @details Starts from the affine map of the bounding box,
   * @f$ \xi_i = 2\,(X_i - \tfrac12(X_{\max,i}+X_{\min,i}))/(X_{\max,i}-X_{\min,i}) @f$, which is exact for an
   * axis-aligned box-shaped cell, and applies Newton's method,
   * @f$ \boldsymbol{\xi} \leftarrow \boldsymbol{\xi} + \boldsymbol{J}^{-1}(\boldsymbol{X} - \sum_A
   * N_A(\boldsymbol{\xi})\boldsymbol{X}_A) @f$, until the residual is below @f$ 10^{-12} @f$ times the diagonal of the
   * bounding box.
   * @param[in] coord Physical coordinates.
   * @return Parametric coordinates @f$ \boldsymbol{\xi} @f$.
   * @throws std::runtime_error if Newton's method does not converge within 20 iterations.
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

  if ( _boundingBoxMatchesGeometryExactly )
    return true;

  XiSized xi;
  try {
    xi = findReferenceCoordinate( Eigen::Map< const XiSized >( coordinates ) );
  }
  catch ( const std::runtime_error& ) {
    return false; // far outside a strongly distorted cell
  }
  for ( auto i = 0; i < nDim; i++ )
    if ( xi( i ) < -1 || xi( i ) >= 1 )
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
  const auto    X = ParentLagrangianGeometryElement::coordinates.reshaped( nDim, nNodes );
  const XiSized h = _boundingBoxMax - _boundingBoxMin;

  // initial guess: the affine map of the bounding box, exact for a box-shaped cell
  XiSized xi = 2 * ( coord - ( _boundingBoxMax + _boundingBoxMin ) / 2 ).cwiseProduct( h.cwiseInverse() );
  XiSized r  = coord - X * N( xi ).transpose();

  for ( int iteration = 0; r.norm() > 1e-12 * h.norm(); iteration++ ) {
    if ( iteration >= 20 )
      throw std::runtime_error( MakeString()
                                << __PRETTY_FUNCTION__ << ": failed to determine inverse map for coordinate "
                                << coord.transpose() );

    const JacobianSized J = ParentLagrangianGeometryElement::Jacobian( ParentLagrangianGeometryElement::dNdXi( xi ) );
    xi += J.partialPivLu().solve( r );
    r = coord - X * N( xi ).transpose();
  }

  return xi;
}
