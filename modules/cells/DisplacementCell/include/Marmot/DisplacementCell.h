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
 * Research Group for Computational Mechanics of Materials
 * Institute of Structural Engineering, BOKU University, Vienna
 *
 * Thomas Mader thomas.mader@boku.ac.at
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
#include "Marmot/DisplacementMaterialPoint.h"
#include "Marmot/MarmotBSplineCellGeometry.h"
#include "Marmot/MarmotCell.h"
#include "Marmot/MarmotCellGeometry.h"
#include "Marmot/MarmotDofLayoutTools.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotLagrangianCellGeometry.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotUtils.h"
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Marmot::Cells {

  /**
   * @class DisplacementCell
   * @brief MPM background cell for finite-strain displacement material points.
   *
   * The cell carries the displacement field (nDim dofs per node) and consumes DisplacementMaterialPoint
   * instances, which drive a MarmotMaterialFiniteStrain. As usual in MPM, the grid dofs are the INCREMENTS of the
   * current step and the material points carry the accumulated state; the momentum balance is assembled with the
   * Kirchhoff stress in the current configuration,
   * @f$ r_{U,Ai} = \sum_p \partial_{x_j} N_A\,\tau_{ij}\,V_p^0 @f$, and its tangent includes the geometric
   * stiffness. It is the displacement-only sibling of GradientEnhancedFiniteStrainCell.
   *
   * @tparam nDim         Spatial dimension (2: plane strain, 3).
   * @tparam nNodes       Number of cell nodes.
   * @tparam CellBase     Base class, MarmotCell.
   * @tparam GeometryCell Geometry policy, Lagrangian or B-spline.
   */
  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  class DisplacementCell : public CellBase, public GeometryCell {

  protected:
    enum BodyLoadTypes {
      BodyForce,
    };

    enum DistributedLoadTypes { Pressure };

    static inline const std::map< std::string, std::pair< int, int > > _fields = {
      { "displacement", { nDim, nNodes } },
    };

    static inline const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };

    static inline const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = {
      { "PRESSURE", Pressure } };

    static constexpr int nDofPerNodeU = nDim; // displacement field U

    static constexpr int bsU            = nNodes * nDofPerNodeU;
    static constexpr int sizeLoadVector = bsU;
    static constexpr int idxU           = 0;

    using MaterialPoint = MaterialPoints::DisplacementMaterialPoint< nDim >;

    using NSized    = typename GeometryCell::NSized;
    using dNdXSized = typename GeometryCell::dNdXSized;
    using XiSized   = typename GeometryCell::XiSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    const int _cellLabel;

    struct MaterialPointLocation {
      MaterialPoint*                         materialPoint;
      XiSized                                xi;
      Fastor::Tensor< double, nNodes >       N;
      Fastor::Tensor< double, nDim, nNodes > dN_dY;
    };

    std::vector< MaterialPointLocation > _materialPointLocations;

  public:
    DisplacementCell( int cellLabel, const GeometryCell& geometry )
      : GeometryCell( geometry ), _cellLabel( cellLabel ){};

    const std::vector< std::vector< std::string > >& getNodeFields() const;

    const std::vector< int >& getDofIndicesPermutationPattern() const;

    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const { return _supportedBodyLoadTypes; }

    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const
    {
      return _supportedDistributedLoadTypes;
    }

    int getNNodes() const { return nNodes; }

    int getNDofPerCell() const { return sizeLoadVector; }

    std::string getCellShape() const { return GeometryCell::getElementShape(); }

    bool isCoordinateInCell( const double* coordinates ) const
    {
      return GeometryCell::isCoordinateInCell( coordinates );
    }

    void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const
    {
      GeometryCell::getBoundingBox( boundingBoxMin, boundingBoxMax );
    }

    void assignMaterialPoints( const std::vector< MarmotMaterialPoint* >& materialPoints );

    void computeMaterialPointKernels( const double* dQ,
                                      double*       fInt,
                                      double*       dfInt_dQ,
                                      double        timeNew,
                                      double        dT ) const;

    void computeLumpedInertia( double* I );

    void computeConsistentInertia( double* I );

    void computeBodyLoad( int type, const double* load, double* fExt, double* dfExt_dQ, double timeNew, double dT )
      const;

    void computeDistributedLoad( int           type,
                                 int           surfaceID,
                                 int           materialPointNumber,
                                 const double* load,
                                 double*       fExt,
                                 double*       dExt_dQ,
                                 double        timeNew,
                                 double        dT ) const;

    void interpolateFieldsToMaterialPoints( const double* dQ ) const;

    void getInterpolationVector( double* vec, const double* coordinates ) const;
  };

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::assignMaterialPoints(
    const std::vector< MarmotMaterialPoint* >& materialPoints )
  {
    _materialPointLocations.clear();

    XiSized coordsMP;

    for ( auto& mp : materialPoints ) {

      auto geMp = dynamic_cast< MaterialPoint* >( mp );
      if ( !geMp )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": material point " << mp->getMaterialPointNumber()
                                     << " is not a DisplacementMaterialPoint" );

      mp->getCoordinatesAtCenter( coordsMP.data() );
      const XiSized xi = GeometryCell::findReferenceCoordinate( coordsMP );

      const auto N_     = GeometryCell::N( xi );
      const auto dN_dY_ = GeometryCell::dNdX( xi );

      _materialPointLocations.push_back( {
        .materialPoint = geMp,
        .xi            = xi,
        .N             = Fastor::Tensor< double, nNodes >( N_.data() ),
        .dN_dY         = Fastor::Tensor< double, nDim, nNodes >( dN_dY_.data(), Fastor::ColumnMajor ),
      } );
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  const std::vector< std::vector< std::string > >& DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::
    getNodeFields() const
  {
    static std::vector< std::vector< std::string > > nodeFields;

    if ( nodeFields.empty() )
      nodeFields = FiniteElement::makeNodeFieldLayout( _fields );

    return nodeFields;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  const std::vector< int >& DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::getDofIndicesPermutationPattern()
    const
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() )
      permutationPattern = FiniteElement::makeBlockedLayoutPermutationPattern( getNodeFields(), _fields );

    return permutationPattern;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::interpolateFieldsToMaterialPoints(
    const double* dQ ) const
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;

    const auto dQU = TensorMap< const double, nNodes, nDim >( dQ );

    for ( auto& mpl : _materialPointLocations ) {

      const auto du    = evaluate( einsum< A, Ai >( mpl.N, dQU ) );
      const auto du_dY = evaluate( einsum< Ai, jA >( dQU, mpl.dN_dY ) );

      mpl.materialPoint->incrementDeformation( du, du_dY );
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::computeMaterialPointKernels( const double* dQ,
                                                                                              double*       fInt_,
                                                                                              double*       dfInt_dQ_,
                                                                                              double        timeNew,
                                                                                              double        dT ) const
  {
    using namespace Fastor;
    using namespace Marmot::FastorIndices;

    Tensor< double, nNodes, nDim >               r_U( 0.0 );
    Tensor< double, nDim, nNodes, nDim, nNodes > k_UU( 0.0 );

    for ( const auto& mpl : _materialPointLocations ) {

      const auto& mp    = mpl.materialPoint;
      const auto& dN_dY = mpl.dN_dY;

      const auto dN_dx = evaluate( einsum< ji, jA >( inv( mp->dx_dY() ), dN_dY ) );

      const double V0 = mp->getVolumeUndeformed();
      const auto&  S  = mp->response.S;

      const auto dS_dqU = evaluate( einsum< ijkl, lB >( mp->tangents.dS_dDeltaF, dN_dY ) );

      r_U += einsum< iA, ij >( dN_dx, S ) * V0;
      k_UU += ( einsum< iA, ijkB, to_jAkB >( dN_dx, dS_dqU ) - einsum< kA, ij, iB, to_jAkB >( dN_dx, S, dN_dx ) ) * V0;
    }

    using namespace Eigen;

    // Due to Fastor Bug #139, we cannot directly write using a TensorMap
    Map< RhsSized >( fInt_ ) += Map< Matrix< double, bsU, 1 > >( r_U.data() );
    Map< KeSizedMatrix >( dfInt_dQ_ ) += Map< Matrix< double, bsU, bsU > >( torowmajor( k_UU ).data() );
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::computeConsistentInertia( double* I )
  {
    Eigen::Map< KeSizedMatrix > M( I );
    M.setZero();

    for ( const auto& mpl : _materialPointLocations ) {
      const double m = mpl.materialPoint->getDensityUndeformed() * mpl.materialPoint->getVolumeUndeformed();
      for ( int A = 0; A < nNodes; A++ )
        for ( int B = 0; B < nNodes; B++ )
          for ( int i = 0; i < nDim; i++ )
            M( idxU + A * nDim + i, idxU + B * nDim + i ) += mpl.N( A ) * mpl.N( B ) * m;
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::computeLumpedInertia( double* I )
  {
    // dynamic storage: a 3D cubic B-spline cell has 256 dofs, too many for a fixed-size matrix on the stack
    Eigen::MatrixXd M( sizeLoadVector, sizeLoadVector );
    computeConsistentInertia( M.data() );
    Eigen::Map< RhsSized > lumped( I );
    lumped = M.rowwise().sum();
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::computeBodyLoad( int           type,
                                                                                  const double* load_,
                                                                                  double*       rhs_,
                                                                                  double*       dRhs_dQ_,
                                                                                  double        timeNew,
                                                                                  double        dT ) const
  {
    switch ( type ) {

    case BodyForce: {

      Fastor::TensorMap< double, nNodes, nDim > r_U( rhs_ );
      Fastor::Tensor< double, nDim >            f( load_ );

      for ( const auto& mpl : _materialPointLocations ) {
        const auto b = Fastor::evaluate( f * mpl.materialPoint->getVolumeUndeformed() );
        r_U -= Fastor::outer( mpl.N, b );
      }
      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid BodyLoadType specified" );
    }
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::computeDistributedLoad( int type,
                                                                                         int surfaceID,
                                                                                         int materialPointNumber,
                                                                                         const double* load_,
                                                                                         double*       fExt_,
                                                                                         double*       dFExt_dQ_,
                                                                                         double        timeNew,
                                                                                         double        dT ) const
  {
    switch ( type ) {

    case Pressure: {

      using namespace Fastor;
      using namespace FastorIndices;

      TensorMap< double, nNodes, nDim > r_U( fExt_ );

      const Tensor< double, nDim > f_0( load_ ); // undeformed load vector p * N * dA_0

      for ( const auto& mpl : _materialPointLocations ) {

        if ( mpl.materialPoint->getMaterialPointNumber() != materialPointNumber )
          continue;

        Tensor< double, nDim, nDim > Eye;
        Eye.eye();

        // Nanson's formula
        const Tensor< double, nDim, nDim > F    = mpl.materialPoint->dx_dY() % mpl.materialPoint->dY_dX();
        const Tensor< double, nDim, nDim > FInv = inverse( F );
        const double                       J    = determinant( F );

        const Tensor< double, nDim > f = J * transpose( FInv ) % f_0;

        const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( FInv, FInv );

        const Tensor< double, nDim, nDim, nDim > df_dF = outer( f, transpose( FInv ) ) +
                                                         J * einsum< IikK, Index< I_ > >( dFInv_dF, f_0 );

        const Tensor< double, nDim, nDim >             F_n        = mpl.materialPoint->dY_dX();
        const Tensor< double, nDim, nDim, nDim, nDim > dF_dDeltaF = einsum< ij, JI, to_iIjJ >( Eye, F_n );

        const Tensor< double, nDim, nDim, nDim >
          df_dDeltaF = einsum< Index< 0, 1, 2 >, Index< 1, 2, 3, 4 > >( df_dF, dF_dDeltaF );
        const Tensor< double, nDim, nDim, nNodes > df_dQU = einsum< Index< 0, 3, 4 >, Index< 4, 5 > >( df_dDeltaF,
                                                                                                       mpl.dN_dY );

        r_U -= outer( mpl.N, f );

        const Tensor< double, nDim, nNodes, nDim, nNodes > dRU_dQU = -einsum< A, jkB, to_jAkB >( mpl.N, df_dQU );

        Eigen::Map< KeSizedMatrix > K( dFExt_dQ_ );
        K.template block< bsU, bsU >( idxU, idxU ) += Eigen::Map< const Eigen::Matrix< double, bsU, bsU > >(
          torowmajor( dRU_dQU ).data() );
      }
      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoadType specified" );
    }
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void DisplacementCell< nDim, nNodes, CellBase, GeometryCell >::getInterpolationVector(
    double*       N,
    const double* coordinates ) const
  {
    const auto refCoord = GeometryCell::findReferenceCoordinate( XiSized( coordinates ) );

    Eigen::Map< Eigen::Matrix< double, nNodes, 1 > > interpolationVector( N );
    interpolationVector = GeometryCell::N( refCoord ).transpose();
  }

  template < int nDim, int nNodes >
  class LagrangianDisplacementCell
    : public DisplacementCell< nDim, nNodes, MarmotCell, MarmotLagrangianCellGeometry< nDim, nNodes > > {

    using Geometry      = MarmotLagrangianCellGeometry< nDim, nNodes >;
    using PhysicsParent = DisplacementCell< nDim, nNodes, MarmotCell, Geometry >;

  public:
    LagrangianDisplacementCell( int cellLabel, const double* nodeCoordinates, int sizeNodeCoordinates )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates ) ){};
  };

  template < int nDim, int nNodes, int order >
  class BSplineDisplacementCell
    : public DisplacementCell< nDim, nNodes, MarmotCell, MarmotBSplineCellGeometry< nDim, order > > {

    using Geometry      = MarmotBSplineCellGeometry< nDim, order >;
    using PhysicsParent = DisplacementCell< nDim, nNodes, MarmotCell, Geometry >;

  public:
    BSplineDisplacementCell( int           cellLabel,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors ) ){};
  };

} // namespace Marmot::Cells
