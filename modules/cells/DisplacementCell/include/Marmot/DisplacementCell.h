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
   * @class Marmot::Cells::DisplacementCell
   * @brief MPM background cell for finite-strain displacement material points.
   *
   * @details The cell carries the displacement field (nDim dofs per node) and consumes DisplacementMaterialPoint
   * instances, which drive a MarmotMaterialFiniteStrain. As usual in MPM, the grid dofs @f$ \Delta q_{Bk} @f$ are the
   * INCREMENTS of the current step, and the material points carry the accumulated state. The cell nodes define the
   * intermediate configuration @f$ \boldsymbol{Y} @f$ (the configuration at the beginning of the increment); the
   * shape functions @f$ N_A @f$ and their gradients @f$ \partial N_A/\partial Y_J @f$ are evaluated once per
   * material point in assignMaterialPoints(). Each material point receives
   * @f[
   *   \Delta\boldsymbol{u} = N_B\,\Delta\boldsymbol{q}_B, \qquad
   *   \Delta F_{iJ} = \delta_{iJ} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial Y_J}.
   * @f]
   * The momentum balance is assembled with the Kirchhoff stress, the spatial gradients
   * @f$ \partial N_A/\partial x_i = \Delta F^{-1}_{Ji}\,\partial N_A/\partial Y_J @f$, and the undeformed volumes
   * @f$ V_p^0 @f$ of the material points,
   * @f[
   *   r_{Aj} = \sum_p \frac{\partial N_A}{\partial x_i}\,\tau_{ij}\,V_p^0,
   * @f]
   * with the tangent (material part and geometric stiffness)
   * @f[
   *   \frac{\partial r_{Aj}}{\partial \Delta q_{Bk}} = \sum_p \left(
   *     \frac{\partial N_A}{\partial x_i}\,\frac{\partial \tau_{ij}}{\partial \Delta F_{kL}}\,
   *     \frac{\partial N_B}{\partial Y_L}
   *     - \frac{\partial N_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i} \right) V_p^0 .
   * @f]
   * It is the displacement-only sibling of GradientEnhancedFiniteStrainCell.
   *
   * @tparam nDim         Spatial dimension (2: plane strain, 3).
   * @tparam nNodes       Number of cell nodes.
   * @tparam CellBase     Base class, MarmotCell.
   * @tparam GeometryCell Geometry policy, Lagrangian or B-spline.
   */
  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  class DisplacementCell : public CellBase, public GeometryCell {

  protected:
    /// Supported body loads.
    enum BodyLoadTypes {
      BodyForce, ///< body force per undeformed volume
    };

    /// Supported distributed loads.
    enum DistributedLoadTypes {
      Pressure ///< surface load at a material point, transformed with Nanson's formula
    };

    /// the node fields: name and (dofs per node, number of nodes)
    static inline const std::map< std::string, std::pair< int, int > > _fields = {
      { "displacement", { nDim, nNodes } },
    };

    /// the names of the supported body loads
    static inline const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };

    /// the names of the supported distributed loads
    static inline const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = {
      { "PRESSURE", Pressure } };

    static constexpr int nDofPerNodeU = nDim;                    ///< dofs per node of the displacement field

    static constexpr int bsU            = nNodes * nDofPerNodeU; ///< size of the displacement block
    static constexpr int sizeLoadVector = bsU;                   ///< number of dofs of the cell
    static constexpr int idxU           = 0;                     ///< first index of the displacement block

    using MaterialPoint = MaterialPoints::DisplacementMaterialPoint< nDim >; ///< the consumed material point type

    using NSized    = typename GeometryCell::NSized;                         ///< shape function vector of the geometry
    using dNdXSized = typename GeometryCell::dNdXSized; ///< shape function gradient matrix of the geometry
    using XiSized   = typename GeometryCell::XiSized;   ///< parametric coordinates

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;              ///< residual vector
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >; ///< stiffness matrix

    const int _cellLabel;                                                          ///< label of the cell

    /**
     * @struct Marmot::Cells::DisplacementCell::MaterialPointLocation
     * @brief A material point hosted by the cell, with its shape functions cached at assignment.
     */
    struct MaterialPointLocation {
      MaterialPoint*                         materialPoint; ///< the material point
      XiSized                                xi;            ///< its parametric coordinates in the cell
      Fastor::Tensor< double, nNodes >       N;             ///< shape functions @f$ N_A @f$
      Fastor::Tensor< double, nDim, nNodes > dN_dY;         ///< gradients @f$ \partial N_A/\partial Y_J @f$
    };

    std::vector< MaterialPointLocation > _materialPointLocations; ///< the currently hosted material points

  public:
    /**
     * @brief Constructs a cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] geometry Geometry of the cell (copied).
     */
    DisplacementCell( int cellLabel, const GeometryCell& geometry )
      : GeometryCell( geometry ), _cellLabel( cellLabel ){};

    /**
     * @brief Node fields of the cell.
     * @return "displacement" at each node.
     */
    const std::vector< std::vector< std::string > >& getNodeFields() const;

    /**
     * @brief Permutation from the node-wise dof ordering to the field-wise ordering.
     * @return The permutation pattern (the identity, as there is a single field).
     */
    const std::vector< int >& getDofIndicesPermutationPattern() const;

    /**
     * @brief Supported body loads.
     * @return "BODYFORCE".
     */
    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const { return _supportedBodyLoadTypes; }

    /**
     * @brief Supported distributed loads.
     * @return "PRESSURE".
     */
    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const
    {
      return _supportedDistributedLoadTypes;
    }

    /**
     * @brief Number of nodes.
     * @return nNodes.
     */
    int getNNodes() const { return nNodes; }

    /**
     * @brief Number of dofs.
     * @return nNodes * nDim.
     */
    int getNDofPerCell() const { return sizeLoadVector; }

    /**
     * @brief Shape of the cell, from the geometry.
     * @return The shape name.
     */
    std::string getCellShape() const { return GeometryCell::getElementShape(); }

    /**
     * @brief Checks whether a point lies in the cell.
     * @param[in] coordinates Coordinates of the point (nDim values).
     * @return True if the point lies in the cell.
     */
    bool isCoordinateInCell( const double* coordinates ) const
    {
      return GeometryCell::isCoordinateInCell( coordinates );
    }

    /**
     * @brief Axis-aligned bounding box of the cell.
     * @param[out] boundingBoxMin Minimum coordinates (nDim values).
     * @param[out] boundingBoxMax Maximum coordinates (nDim values).
     */
    void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const
    {
      GeometryCell::getBoundingBox( boundingBoxMin, boundingBoxMax );
    }

    /**
     * @brief Assigns the material points currently located in the cell and caches, at their positions
     * @f$ \boldsymbol{Y} @f$ (see DisplacementMaterialPoint::getCoordinatesAtCenter()), the shape functions and their
     * gradients.
     * @param[in] materialPoints The material points.
     * @throws std::invalid_argument if a material point is not a DisplacementMaterialPoint of dimension nDim.
     */
    void assignMaterialPoints( const std::vector< MarmotMaterialPoint* >& materialPoints );

    /**
     * @brief Assembles the internal force vector and its tangent of the hosted material points (see the class
     * description). The material points must have been updated by interpolateFieldsToMaterialPoints() and computed
     * before.
     * @param[in] dQ Nodal displacement increments (not used; the kinematics are taken from the material points).
     * @param[in,out] fInt Internal force vector, the contribution is added.
     * @param[in,out] dfInt_dQ Tangent @f$ \partial r/\partial \Delta q @f$, the contribution is added.
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     */
    void computeMaterialPointKernels( const double* dQ,
                                      double*       fInt,
                                      double*       dfInt_dQ,
                                      double        timeNew,
                                      double        dT ) const;

    /**
     * @brief Lumped mass vector, the row sums of the consistent mass matrix.
     * @param[out] I Lumped mass vector (sizeLoadVector values), overwritten.
     */
    void computeLumpedInertia( double* I );

    /**
     * @brief Consistent mass matrix @f$ M_{AiBi} = \sum_p N_A\,N_B\,\rho_0\,V_p^0 @f$.
     * @param[out] I Mass matrix (sizeLoadVector x sizeLoadVector), overwritten.
     */
    void computeConsistentInertia( double* I );

    /**
     * @brief Body load: @f$ r_{Aj} \mathrel{-}= \sum_p N_A\,f_j\,V_p^0 @f$ with a body force @f$ \boldsymbol{f} @f$
     * per undeformed volume; the load does not contribute to the tangent.
     * @param[in] type The body load type (BodyForce).
     * @param[in] load The body force (nDim values).
     * @param[in,out] fExt Load vector, the contribution is added.
     * @param[in,out] dfExt_dQ Tangent (not modified).
     * @param[in] timeNew Time at the end of the increment.
     * @param[in] dT Time increment.
     * @throws std::invalid_argument for an unsupported type.
     */
    void computeBodyLoad( int type, const double* load, double* fExt, double* dfExt_dQ, double timeNew, double dT )
      const;

    /**
     * @brief Distributed load at a single material point of the cell.
     *
     * @details The load vector @f$ \boldsymbol{f}_0 @f$ (e.g. @f$ p\,\boldsymbol{N}\,dA_0 @f$ in the undeformed
     * configuration) is transformed with Nanson's formula and the total deformation gradient,
     * @f$ \boldsymbol{f} = J\,\boldsymbol{F}^{-\mathsf T}\boldsymbol{f}_0 @f$, and assembled as
     * @f$ r_{Aj} \mathrel{-}= N_A\,f_j @f$, with the corresponding load stiffness.
     * @param[in] type The distributed load type (Pressure).
     * @param[in] surfaceID Surface ID (not used).
     * @param[in] materialPointNumber Label of the material point the load acts on; other material points are skipped.
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

    /**
     * @brief Interpolates the nodal increments to the hosted material points, see the class description; calls
     * DisplacementMaterialPoint::incrementDeformation() (which accumulates, so the host resets the material points
     * with DisplacementMaterialPoint::prepareYourself() first).
     * @param[in] dQ Nodal displacement increments (node-wise, nDim values per node).
     */
    void interpolateFieldsToMaterialPoints( const double* dQ ) const;

    /**
     * @brief Shape functions of the cell at an arbitrary point.
     * @param[out] vec Shape functions (nNodes values).
     * @param[in] coordinates Coordinates of the point (nDim values).
     */
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

    Tensor< double, nNodes, nDim > r_U( 0.0 );

    // the tangent is assembled directly into the output: a fixed-size nDim*nNodes x nDim*nNodes tensor and its einsum
    // temporaries (~300 kB each for a 64-node hexahedron) overflow the 1 MB default stack of Windows
    Eigen::Map< Eigen::MatrixXd > K( dfInt_dQ_, sizeLoadVector, sizeLoadVector );

    for ( const auto& mpl : _materialPointLocations ) {

      const auto& mp    = mpl.materialPoint;
      const auto& dN_dY = mpl.dN_dY;

      const auto dN_dx = evaluate( einsum< ji, jA >( inv( mp->dx_dY() ), dN_dY ) );

      const double V0 = mp->getVolumeUndeformed();
      const auto&  S  = mp->response.S;

      const auto dS_dqU = evaluate( einsum< ijkl, lB >( mp->tangents.dS_dDeltaF, dN_dY ) );

      r_U += einsum< iA, ij >( dN_dx, S ) * V0;

      // K_(Aj)(Bk) = ( dN_A/dx_i dS_ij/dq_Bk - dN_A/dx_k S_ij dN_B/dx_i ) V0
      const auto SdN_dx = evaluate( einsum< ij, iB >( S, dN_dx ) ); // S_ij dN_B/dx_i
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

    using namespace Eigen;

    // Due to Fastor Bug #139, we cannot directly write using a TensorMap
    Map< RhsSized >( fInt_ ) += Map< Matrix< double, bsU, 1 > >( r_U.data() );
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

        // dR_(Aj)/dQ_(Bk) = -N_A df_j/dQ_Bk, assembled directly: a fixed-size tensor of the full block (~300 kB for a
        // 64-node hexahedron) and its temporaries overflow the 1 MB default stack of Windows
        Eigen::Map< Eigen::MatrixXd > K( dFExt_dQ_, sizeLoadVector, sizeLoadVector );
        for ( int A = 0; A < nNodes; A++ )
          for ( int j = 0; j < nDim; j++ )
            for ( int B = 0; B < nNodes; B++ )
              for ( int k = 0; k < nDim; k++ )
                K( idxU + A * nDim + j, idxU + B * nDim + k ) -= mpl.N( A ) * df_dQU( j, k, B );
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

  /**
   * @class Marmot::Cells::LagrangianDisplacementCell
   * @brief DisplacementCell with a Lagrangian geometry (registered as "Displacement/Quad4" and "Displacement/Hexa8").
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of nodes (4 in 2D, 8 in 3D).
   */
  template < int nDim, int nNodes >
  class LagrangianDisplacementCell
    : public DisplacementCell< nDim, nNodes, MarmotCell, MarmotLagrangianCellGeometry< nDim, nNodes > > {

    using Geometry      = MarmotLagrangianCellGeometry< nDim, nNodes >;           ///< the geometry policy
    using PhysicsParent = DisplacementCell< nDim, nNodes, MarmotCell, Geometry >; ///< the physics base class

  public:
    /**
     * @brief Constructs a Lagrangian cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] nodeCoordinates Node coordinates (nDim values per node).
     * @param[in] sizeNodeCoordinates Number of coordinates (not used).
     */
    LagrangianDisplacementCell( int cellLabel, const double* nodeCoordinates, int sizeNodeCoordinates )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates ) ){};
  };

  /**
   * @class Marmot::Cells::BSplineDisplacementCell
   * @brief DisplacementCell with a B-spline geometry (registered as "Displacement/BSpline/<order>" in 2D and
   * "Displacement/BSpline/3D/<order>" in 3D).
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of control points, @f$ (order+1)^{nDim} @f$.
   * @tparam order  Polynomial order of the B-splines (1, 2 or 3).
   */
  template < int nDim, int nNodes, int order >
  class BSplineDisplacementCell
    : public DisplacementCell< nDim, nNodes, MarmotCell, MarmotBSplineCellGeometry< nDim, order > > {

    using Geometry      = MarmotBSplineCellGeometry< nDim, order >;               ///< the geometry policy
    using PhysicsParent = DisplacementCell< nDim, nNodes, MarmotCell, Geometry >; ///< the physics base class

  public:
    /**
     * @brief Constructs a B-spline cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] nodeCoordinates Control point coordinates (nDim values per control point).
     * @param[in] sizeNodeCoordinates Number of coordinates.
     * @param[in] knotVectors Knot vectors of the cell.
     * @param[in] sizeKnotVectors Number of knot values.
     */
    BSplineDisplacementCell( int           cellLabel,
                             const double* nodeCoordinates,
                             int           sizeNodeCoordinates,
                             const double* knotVectors,
                             int           sizeKnotVectors )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors ) ){};
  };

} // namespace Marmot::Cells
