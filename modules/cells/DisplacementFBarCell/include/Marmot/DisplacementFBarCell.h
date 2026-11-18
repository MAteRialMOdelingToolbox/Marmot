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
#include "Marmot/DisplacementCell.h"
#include <cmath>
#include <type_traits>
#include <vector>

namespace Marmot::Discretizations::Cells {

  /**
   * @class Marmot::Discretizations::Cells::DisplacementFBarCell
   * @brief MPM background cell for finite-strain displacement material points with an F-bar projection of the
   * volumetric deformation.
   *
   * @details The cell is a DisplacementCell (same dofs, loads, mass and geometry) whose material points see a
   * modified incremental deformation gradient. The volume change of the increment is taken at the geometric center
   * @f$ \bar{\boldsymbol{\xi}} = \frac{1}{n_p}\sum_p \boldsymbol{\xi}_p @f$ of the hosted material points (the mean of
   * their parametric coordinates, which for an affine cell is the image of their centroid),
   * @f[
   *   \Delta\boldsymbol{F}_c = \boldsymbol{I} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial \boldsymbol{Y}}
   *     \Big|_{\bar{\boldsymbol{\xi}}}, \qquad
   *   \Delta\bar{\boldsymbol{F}}_p = s_p\,\Delta\boldsymbol{F}_p, \qquad
   *   s_p = \left( \frac{\det \Delta\boldsymbol{F}_c}{\det \Delta\boldsymbol{F}_p} \right)^{1/n_\mathrm{dim}},
   * @f]
   * and every material point receives @f$ \Delta\bar{\boldsymbol{F}}_p @f$ instead of
   * @f$ \Delta\boldsymbol{F}_p @f$ (through DisplacementMaterialPoint::incrementDeformation(); its accumulated
   * deformation gradient is hence @f$ \bar{\boldsymbol{F}} = \Delta\bar{\boldsymbol{F}}\,\bar{\boldsymbol{F}}_n @f$).
   * The projection is incremental because the cell has no history of its own: the material points carry the
   * accumulated state. In plane strain the exponent is @f$ 1/2 @f$ and the out-of-plane stretch stays 1.
   *
   * The internal forces are those of the DisplacementCell with the Kirchhoff stress
   * @f$ \boldsymbol{\tau}(\bar{\boldsymbol{F}}) @f$ and the spatial gradients of the *unmodified* kinematics,
   * @f$ \partial N_A/\partial x_i = \Delta F^{-1}_{p,Ji}\,\partial N_A/\partial Y_J @f$. The tangent is exact: with
   * @f$ \partial \ln\det\Delta\boldsymbol{F}_p / \partial \Delta q_{Bm} = \partial N_B/\partial x_m @f$ and the
   * analogous center gradient @f$ \partial N^c_B/\partial x_m = \Delta F^{-1}_{c,Lm}\,\partial N_B/\partial Y_L
   * |_{\bar{\boldsymbol{\xi}}} @f$,
   * @f[
   *   \frac{\partial \Delta\bar F_{kL}}{\partial \Delta q_{Bm}}
   *     = s_p\,\delta_{km}\,\frac{\partial N_B}{\partial Y_L}
   *     + \frac{1}{n_\mathrm{dim}}\,\Delta\bar F_{kL}
   *       \left( \frac{\partial N^c_B}{\partial x_m} - \frac{\partial N_B}{\partial x_m} \right),
   * @f]
   * which replaces @f$ \delta_{km}\,\partial N_B/\partial Y_L @f$ in the material part of the DisplacementCell tangent;
   * the geometric part is unchanged. The follower pressure load is transformed with
   * @f$ \bar{\boldsymbol{F}} @f$, the deformation gradient of the material point, with the same chain rule.
   *
   * With a single material point in the cell, the center is the point itself, @f$ s_p \equiv 1 @f$, and the cell
   * coincides with the DisplacementCell (residual and tangent). For a homogeneous increment the residual coincides
   * with that of the DisplacementCell.
   *
   * @tparam nDim         Spatial dimension (2: plane strain, 3).
   * @tparam nNodes       Number of cell nodes.
   * @tparam CellBase     Base class, Marmot::Discretizations::Cells::Cell.
   * @tparam GeometryCell Geometry policy, Lagrangian or B-spline.
   */
  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  class DisplacementFBarCell : public DisplacementCell< nDim, nNodes, CellBase, GeometryCell > {

  protected:
    using Parent   = DisplacementCell< nDim, nNodes, CellBase, GeometryCell >; ///< the displacement cell
    using XiSized  = typename Parent::XiSized;                                 ///< parametric coordinates
    using TensorDD = Fastor::Tensor< double, nDim, nDim >;                     ///< second-order tensor of size nDim
    using dFdQ = Fastor::Tensor< double, nDim, nDim, nDim, nNodes >; ///< derivative of a deformation gradient w.r.t.
                                                                     ///< the nodal displacement increments

    /**
     * @struct Marmot::Discretizations::Cells::DisplacementFBarCell::Center
     * @brief The geometric center of the hosted material points.
     */
    struct Center {
      XiSized                                xi;    ///< its parametric coordinates
      Fastor::Tensor< double, nDim, nNodes > dN_dY; ///< shape function gradients @f$ \partial N_B/\partial Y_J @f$
    };

    Center _center;                        ///< the geometric center of the hosted material points

    mutable TensorDD                _dF_c; ///< @f$ \Delta\boldsymbol{F}_c @f$ of the last interpolation
    mutable std::vector< TensorDD > _dF_p; ///< unmodified @f$ \Delta\boldsymbol{F}_p @f$ of the last interpolation

    /**
     * @brief The derivative of the modified increment of a material point w.r.t. the nodal displacement increments,
     * see the class description.
     * @param[in] p Index of the material point in the cell.
     * @return @f$ \partial \Delta\bar F_{kL}/\partial \Delta q_{Bm} @f$, indexed (k, L, m, B).
     */
    dFdQ dDeltaFBar_dQ( size_t p ) const;

  public:
    using Parent::Parent;

    /**
     * @brief Assigns the material points (see DisplacementCell::assignMaterialPoints()) and evaluates the shape
     * function gradients at their geometric center.
     * @param[in] materialPoints The material points.
     */
    void assignMaterialPoints(
      const std::vector< Marmot::Discretizations::MaterialPoints::MaterialPoint* >& materialPoints );

    /**
     * @brief Evaluates the volume change of the cell at a given point instead of the geometric center of its material
     * points, until the next assignMaterialPoints(); see Cell::setProjectionPoint(). Outside the cell, the cell's
     * basis is extended (aggregated MPM: all cells of an aggregate evaluate @f$ \Delta\boldsymbol{F}_c @f$ at the same
     * point of the field they share).
     * @param[in] coordinates The point, in physical coordinates.
     * @return true for Lagrangian geometries; false for B-spline geometries, whose gradients are not extended beyond
     * the cell (the center of the material points is kept).
     */
    bool setProjectionPoint( const double* coordinates ) override;

    /**
     * @brief Interpolates the nodal increments to the hosted material points; they receive the displacement increment
     * and the modified incremental deformation gradient @f$ \Delta\bar{\boldsymbol{F}}_p @f$ (see the class
     * description).
     * @param[in] dQ Nodal displacement increments (node-wise, nDim values per node).
     */
    void interpolateFieldsToMaterialPoints( const double* dQ ) const;

    /**
     * @brief Internal force vector and its exact tangent, see the class description. The material points must have
     * been updated by interpolateFieldsToMaterialPoints() and computed before.
     * @param[in] dQ Nodal displacement increments (not used).
     * @param[in,out] fInt Internal force vector, the contribution is added.
     * @param[in,out] dfInt_dQ Tangent, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    void computeMaterialPointKernels( const double* dQ,
                                      double*       fInt,
                                      double*       dfInt_dQ,
                                      double        timeNew,
                                      double        dT ) const;

    /**
     * @brief Follower pressure load at a single material point, see DisplacementCell::computeDistributedLoad(), with
     * the modified deformation gradient and its tangent.
     * @param[in] type The distributed load type (Pressure).
     * @param[in] surfaceID Surface ID (not used).
     * @param[in] materialPointNumber Label of the material point the load acts on.
     * @param[in] load The load vector @f$ \boldsymbol{f}_0 @f$ (nDim values).
     * @param[in,out] fExt Load vector, the contribution is added.
     * @param[in,out] dExt_dQ Tangent, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::invalid_argument for an unsupported type.
     */
    void computeDistributedLoad( int           type,
                                 int           surfaceID,
                                 int           materialPointNumber,
                                 const double* load,
                                 double*       fExt,
                                 double*       dExt_dQ,
                                 double        timeNew,
                                 double        dT ) const;
  };

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::assignMaterialPoints(
    const std::vector< Marmot::Discretizations::MaterialPoints::MaterialPoint* >& materialPoints )
  {
    Parent::assignMaterialPoints( materialPoints );

    const auto& locations = this->_materialPointLocations;

    _center.xi = XiSized::Zero();
    for ( const auto& mpl : locations )
      _center.xi += mpl.xi;
    if ( !locations.empty() )
      _center.xi /= static_cast< double >( locations.size() );

    const auto dN_dY_ = GeometryCell::dNdX( _center.xi );
    _center.dN_dY     = Fastor::Tensor< double, nDim, nNodes >( dN_dY_.data(), Fastor::ColumnMajor );

    _dF_c.eye();
    _dF_p.assign( locations.size(), _dF_c );
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  bool DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::setProjectionPoint( const double* coordinates )
  {
    if constexpr ( std::is_same_v< GeometryCell, LagrangianCellGeometry< nDim, nNodes > > ) {
      _center.xi        = GeometryCell::findReferenceCoordinate( Eigen::Map< const XiSized >( coordinates ) );
      const auto dN_dY_ = GeometryCell::dNdX( _center.xi );
      _center.dN_dY     = Fastor::Tensor< double, nDim, nNodes >( dN_dY_.data(), Fastor::ColumnMajor );
      return true;
    }
    else
      return false;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::interpolateFieldsToMaterialPoints(
    const double* dQ ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    const auto dQU = TensorMap< const double, nNodes, nDim >( dQ );

    TensorDD Eye;
    Eye.eye();

    _dF_c               = Eye + einsum< Ai, jA >( dQU, _center.dN_dY );
    const double detF_c = determinant( _dF_c );

    const auto& locations = this->_materialPointLocations;
    for ( size_t p = 0; p < locations.size(); p++ ) {
      const auto& mpl = locations[p];

      const auto du = evaluate( einsum< A, Ai >( mpl.N, dQU ) );
      _dF_p[p]      = Eye + einsum< Ai, jA >( dQU, mpl.dN_dY );

      const double   s         = std::pow( detF_c / determinant( _dF_p[p] ), 1. / nDim );
      const TensorDD du_dY_bar = s * _dF_p[p] - Eye;

      mpl.materialPoint->incrementDeformation( du, du_dY_bar );
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  typename DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::dFdQ DisplacementFBarCell<
    nDim,
    nNodes,
    CellBase,
    GeometryCell >::dDeltaFBar_dQ( size_t p ) const
  {
    using namespace Marmot::Math::Tensors::FastorIndices;
    using namespace Fastor;

    const auto& dN_dY = this->_materialPointLocations[p].dN_dY;
    const auto& dF_p  = _dF_p[p];

    const double   s     = std::pow( determinant( _dF_c ) / determinant( dF_p ), 1. / nDim );
    const TensorDD dFbar = s * dF_p;

    // d ln det(dF)/d dq_Bm, at the material point and at the center
    const Tensor< double, nDim, nNodes > dN_dx  = einsum< ji, jA >( inverse( dF_p ), dN_dY );
    const Tensor< double, nDim, nNodes > dNc_dx = einsum< ji, jA >( inverse( _dF_c ), _center.dN_dY );

    dFdQ result;
    for ( int k = 0; k < nDim; k++ )
      for ( int L = 0; L < nDim; L++ )
        for ( int m = 0; m < nDim; m++ )
          for ( int B = 0; B < nNodes; B++ )
            result( k, L, m, B ) = ( k == m ? s * dN_dY( L, B ) : 0.0 ) +
                                   dFbar( k, L ) / nDim * ( dNc_dx( m, B ) - dN_dx( m, B ) );
    return result;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::computeMaterialPointKernels( const double* dQ,
                                                                                                  double*       fInt_,
                                                                                                  double* dfInt_dQ_,
                                                                                                  double  timeNew,
                                                                                                  double  dT ) const
  {
    using namespace Fastor;
    using namespace Marmot::Math::Tensors::FastorIndices;

    Tensor< double, nNodes, nDim > r_U( 0.0 );

    // assembled directly into the output, see DisplacementCell::computeMaterialPointKernels()
    Eigen::Map< Eigen::MatrixXd > K( dfInt_dQ_, Parent::sizeLoadVector, Parent::sizeLoadVector );

    const auto& locations = this->_materialPointLocations;
    for ( size_t p = 0; p < locations.size(); p++ ) {

      const auto& mp = locations[p].materialPoint;

      // the spatial gradients of the unmodified kinematics
      const auto dN_dx = evaluate( einsum< ji, jA >( inverse( _dF_p[p] ), locations[p].dN_dY ) );

      const double V0 = mp->getVolumeUndeformed();
      const auto&  S  = mp->response.S;
      const auto&  C  = mp->tangents.dS_dDeltaF; // w.r.t. the modified increment

      const dFdQ dFbar_dQ = dDeltaFBar_dQ( p );

      // dS_ij/dq_Bm
      Tensor< double, nDim, nDim, nDim, nNodes > dS_dqU( 0.0 );
      for ( int i = 0; i < nDim; i++ )
        for ( int j = 0; j < nDim; j++ )
          for ( int k = 0; k < nDim; k++ )
            for ( int L = 0; L < nDim; L++ )
              for ( int m = 0; m < nDim; m++ )
                for ( int B = 0; B < nNodes; B++ )
                  dS_dqU( i, j, m, B ) += C( i, j, k, L ) * dFbar_dQ( k, L, m, B );

      r_U += einsum< iA, ij >( dN_dx, S ) * V0;

      // K_(Aj)(Bk) = ( dN_A/dx_i dS_ij/dq_Bk - dN_A/dx_k S_ij dN_B/dx_i ) V0
      const auto SdN_dx = evaluate( einsum< ij, iB >( S, dN_dx ) );
      for ( int A = 0; A < nNodes; A++ )
        for ( int j = 0; j < nDim; j++ )
          for ( int B = 0; B < nNodes; B++ )
            for ( int k = 0; k < nDim; k++ ) {
              double kAjBk = -dN_dx( k, A ) * SdN_dx( j, B );
              for ( int i = 0; i < nDim; i++ )
                kAjBk += dN_dx( i, A ) * dS_dqU( i, j, k, B );
              K( A * nDim + j, B * nDim + k ) += kAjBk * V0;
            }
    }

    Eigen::Map< typename Parent::RhsSized >( fInt_ ) += Eigen::Map< Eigen::Matrix< double, Parent::bsU, 1 > >(
      r_U.data() );
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementFBarCell< nDim, nNodes, CellBase, GeometryCell >::computeDistributedLoad( int type,
                                                                                             int surfaceID,
                                                                                             int materialPointNumber,
                                                                                             const double* load_,
                                                                                             double*       fExt_,
                                                                                             double*       dFExt_dQ_,
                                                                                             double        timeNew,
                                                                                             double        dT ) const
  {
    switch ( type ) {

    case Parent::Pressure: {

      using namespace Fastor;
      using namespace Marmot::Math::Tensors::FastorIndices;

      TensorMap< double, nNodes, nDim > r_U( fExt_ );
      Eigen::Map< Eigen::MatrixXd >     K( dFExt_dQ_, Parent::sizeLoadVector, Parent::sizeLoadVector );

      const Tensor< double, nDim > f_0( load_ ); // undeformed load vector p * N * dA_0

      const auto& locations = this->_materialPointLocations;
      for ( size_t p = 0; p < locations.size(); p++ ) {
        const auto& mpl = locations[p];

        if ( mpl.materialPoint->getMaterialPointNumber() != materialPointNumber )
          continue;

        // Nanson's formula with the deformation gradient of the material point, dFbar F_n
        const TensorDD F_n  = mpl.materialPoint->dY_dX();
        const TensorDD F    = mpl.materialPoint->dx_dY() % F_n;
        const TensorDD FInv = inverse( F );
        const double   detF = determinant( F );

        const Tensor< double, nDim > f = detF * transpose( FInv ) % f_0;

        // df_j/dF_kQ = f_j FInv_Qk - J FInv_Qj FInv_Rk f0_R, and df_j/d dFbar_kL = df_j/dF_kQ F_n,LQ
        Tensor< double, nDim, nDim, nDim > df_dDFbar( 0.0 );
        for ( int j = 0; j < nDim; j++ )
          for ( int k = 0; k < nDim; k++ )
            for ( int L = 0; L < nDim; L++ )
              for ( int Q = 0; Q < nDim; Q++ ) {
                double df_dF = f( j ) * FInv( Q, k );
                for ( int R = 0; R < nDim; R++ )
                  df_dF -= detF * FInv( Q, j ) * FInv( R, k ) * f_0( R );
                df_dDFbar( j, k, L ) += df_dF * F_n( L, Q );
              }

        const dFdQ dFbar_dQ = dDeltaFBar_dQ( p );

        r_U -= outer( mpl.N, f );

        for ( int A = 0; A < nNodes; A++ )
          for ( int j = 0; j < nDim; j++ )
            for ( int B = 0; B < nNodes; B++ )
              for ( int m = 0; m < nDim; m++ ) {
                double df_dq = 0;
                for ( int k = 0; k < nDim; k++ )
                  for ( int L = 0; L < nDim; L++ )
                    df_dq += df_dDFbar( j, k, L ) * dFbar_dQ( k, L, m, B );
                K( Parent::idxU + A * nDim + j, Parent::idxU + B * nDim + m ) -= mpl.N( A ) * df_dq;
              }
      }
      break;
    }
    default: {
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid DistributedLoadType specified" );
    }
    }
  }

  /**
   * @class Marmot::Discretizations::Cells::LagrangianDisplacementFBarCell
   * @brief DisplacementFBarCell with a Lagrangian geometry (registered as "DisplacementFBar/Quad4" and
   * "DisplacementFBar/Hexa8").
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of nodes (4 in 2D, 8 in 3D).
   */
  template < int nDim, int nNodes >
  class LagrangianDisplacementFBarCell
    : public DisplacementFBarCell< nDim,
                                   nNodes,
                                   Marmot::Discretizations::Cells::Cell,
                                   Marmot::Discretizations::Cells::LagrangianCellGeometry< nDim, nNodes > > {

    using Geometry = Marmot::Discretizations::Cells::LagrangianCellGeometry< nDim, nNodes >; ///< the geometry policy
    using PhysicsParent = DisplacementFBarCell< nDim,
                                                nNodes,
                                                Marmot::Discretizations::Cells::Cell,
                                                Geometry >; ///< the physics base class

  public:
    /**
     * @brief Constructs a Lagrangian cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] nodeCoordinates Node coordinates (nDim values per node).
     * @param[in] sizeNodeCoordinates Number of coordinates (not used).
     */
    LagrangianDisplacementFBarCell( int cellLabel, const double* nodeCoordinates, int sizeNodeCoordinates )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates ) ){};
  };

  /**
   * @class Marmot::Discretizations::Cells::BSplineDisplacementFBarCell
   * @brief DisplacementFBarCell with a B-spline geometry (registered as "DisplacementFBar/BSpline/<order>" in 2D and
   * "DisplacementFBar/BSpline/3D/<order>" in 3D).
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of control points, @f$ (order+1)^{nDim} @f$.
   * @tparam order  Polynomial order of the B-splines (1, 2 or 3).
   */
  template < int nDim, int nNodes, int order >
  class BSplineDisplacementFBarCell
    : public DisplacementFBarCell< nDim,
                                   nNodes,
                                   Marmot::Discretizations::Cells::Cell,
                                   Marmot::Discretizations::Cells::BSplineCellGeometry< nDim, order > > {

    using Geometry      = Marmot::Discretizations::Cells::BSplineCellGeometry< nDim, order >; ///< the geometry policy
    using PhysicsParent = DisplacementFBarCell< nDim,
                                                nNodes,
                                                Marmot::Discretizations::Cells::Cell,
                                                Geometry >; ///< the physics base class

  public:
    /**
     * @brief Constructs a B-spline cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] nodeCoordinates Control point coordinates (nDim values per control point).
     * @param[in] sizeNodeCoordinates Number of coordinates.
     * @param[in] knotVectors Knot vectors of the cell.
     * @param[in] sizeKnotVectors Number of knot values.
     */
    BSplineDisplacementFBarCell( int           cellLabel,
                                 const double* nodeCoordinates,
                                 int           sizeNodeCoordinates,
                                 const double* knotVectors,
                                 int           sizeKnotVectors )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors ) ){};
  };

} // namespace Marmot::Discretizations::Cells
