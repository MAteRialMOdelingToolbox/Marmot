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
#include "Marmot/GradientEnhancedFiniteStrainMaterialPoint.h"
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
   * @class Marmot::Cells::GradientEnhancedFiniteStrainCell
   * @brief MPM background cell for gradient-enhanced (implicit-gradient) finite-strain materials.
   *
   * The cell carries the displacement field @f$ \boldsymbol{u} @f$ (nDim dofs per node) and one scalar nonlocal
   * field @f$ \bar{N} @f$ (1 dof per node, field name `nonlocal damage`), and consumes
   * Marmot::MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint instances, which in turn drive a
   * MarmotMaterialGradientEnhancedFiniteStrain. It is the MPM counterpart of the finite element
   * GradientEnhancedFiniteStrainDisplacementElement.
   *
   * As usual in MPM, the grid dofs are the INCREMENTS of the current step, @f$ \Delta\boldsymbol{q} @f$; the
   * material points carry the accumulated state. The cell shape functions @f$ N_A @f$ and their gradients
   * @f$ \nabla_Y N_A @f$ are evaluated once, at the position of each material point in the intermediate reference
   * configuration @f$ \boldsymbol{Y} @f$ (the last accepted configuration), when the material points are assigned.
   * The spatial and the reference gradients follow from
   * @f[
   *   \frac{\partial N_A}{\partial x_i} = \frac{\partial N_A}{\partial Y_j}\,\Delta F^{-1}_{ji}, \qquad
   *   \frac{\partial N_A}{\partial X_i} = \frac{\partial N_A}{\partial Y_j}\,F_{n,ji} .
   * @f]
   * The momentum balance is integrated with the Kirchhoff stress over the undeformed volume @f$ V_p^0 @f$ of each
   * material point @f$ p @f$,
   * @f[
   *   r_{U,Aj} = \sum_p \frac{\partial N_A}{\partial x_i}\,\tau_{ij}\,V_p^0 ,
   * @f]
   * and the nonlocal balance @f$ \bar{N} - c\,\nabla_X^2\bar{N} = L @f$ is assembled in increment form in the
   * undeformed configuration,
   * @f[
   *   r_{N,A} = \sum_p \Bigl( N_A\,\Delta\bar{N}_p + c\,\nabla_X N_A \cdot \nabla_X \Delta\bar{N}_p - N_A\,\Delta L_p
   *   \Bigr) V_p^0 ,
   * @f]
   * with @f$ \Delta\bar{N}_p = N_B\,\Delta\bar{N}_B @f$ interpolated from the grid increments, @f$ \Delta L_p @f$ the
   * change of the local driving force reported by the material point, and @f$ c = R^2 @f$ from the material. The
   * reference gradient @f$ \nabla_X @f$ does not depend on the current increment.
   *
   * The increment form is used because the grid is reset every increment, so its nodal values carry only the
   * increment of the step; the total form of the element of the finite element method,
   * @f$ N_A(\bar{N} - L) + c\,\nabla_X N_A\cdot\nabla_X\bar{N} @f$, would need @f$ \nabla_X\bar{N} @f$ as an
   * accumulated state of the material point. It thus assumes that the equation of the previous increment holds with
   * the current grid and point positions; storing @f$ \nabla_X\bar{N} @f$ instead was tested and changed the
   * nonlocal field by less than the step dependence of the kinematics, without bringing it closer to the result of
   * a single increment.
   *
   * The consistent tangent consists of
   * @f[
   *   K^{UU}_{jAkB} = \sum_p \Bigl( \frac{\partial N_A}{\partial x_i}\,
   *     \frac{\partial\tau_{ij}}{\partial\Delta F_{kL}}\,\frac{\partial N_B}{\partial Y_L}
   *     - \frac{\partial N_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i} \Bigr) V_p^0 ,\qquad
   *   K^{UN}_{jAB} = \sum_p \frac{\partial N_A}{\partial x_i}\,\frac{\partial\tau_{ij}}{\partial\bar{N}}\,N_B\,V_p^0 ,
   * @f]
   * @f[
   *   K^{NU}_{AkB} = -\sum_p N_A\,\frac{\partial L}{\partial\Delta F_{kL}}\,\frac{\partial N_B}{\partial Y_L}\,V_p^0 ,
   *   \qquad
   *   K^{NN}_{AB} = \sum_p \Bigl( N_A N_B \bigl( 1 - \frac{\partial L}{\partial\bar{N}} \bigr)
   *     + c\,\nabla_X N_A\cdot\nabla_X N_B \Bigr) V_p^0 .
   * @f]
   * The interaction @f$ c @f$ is treated as constant in the tangent.
   *
   * Dofs are ordered field by field internally (all displacements, then all nonlocal dofs); see
   * getDofIndicesPermutationPattern() for the mapping to the node-by-node layout of the host.
   *
   * @tparam nDim         Spatial dimension (2: plane strain, 3).
   * @tparam nNodes       Number of cell nodes.
   * @tparam CellBase     Base class, MarmotCell.
   * @tparam GeometryCell Geometry policy, Lagrangian or B-spline.
   */
  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  class GradientEnhancedFiniteStrainCell : public CellBase, public GeometryCell {

  protected:
    /// @brief Supported body loads.
    enum BodyLoadTypes {
      BodyForce, ///< body force per unit undeformed volume (`BODYFORCE`)
    };

    /// @brief Supported distributed loads.
    enum DistributedLoadTypes {
      Pressure ///< follower pressure acting on a material point (`PRESSURE`)
    };

    /// the nodal fields with their number of components and nodes
    static inline const std::map< std::string, std::pair< int, int > > _fields = {
      { "displacement", { nDim, nNodes } },
      { "nonlocal damage", { 1, nNodes } },
    };

    /// names of the supported body loads
    static inline const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };

    /// names of the supported distributed loads
    static inline const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = {
      { "PRESSURE", Pressure } };

    static constexpr int nDofPerNodeU = nDim; ///< dofs per node of the displacement field U
    static constexpr int nDofPerNodeN = 1;    ///< dofs per node of the nonlocal field N

    // block sizes
    static constexpr int bsU = nNodes * nDofPerNodeU; ///< size of the displacement block
    static constexpr int bsN = nNodes * nDofPerNodeN; ///< size of the nonlocal block

    static constexpr int sizeLoadVector = bsU + bsN;  ///< number of dofs of the cell

    static constexpr int idxU = 0;                    ///< first index of the displacement block
    static constexpr int idxN = idxU + bsU;           ///< first index of the nonlocal block

    using MaterialPoint = MaterialPoints::GradientEnhancedFiniteStrainMaterialPoint< nDim >; ///< consumed point type

    using NSized    = typename GeometryCell::NSized;                  ///< shape function vector of the geometry
    using dNdXSized = typename GeometryCell::dNdXSized;               ///< shape function gradients of the geometry
    using XiSized   = typename GeometryCell::XiSized;                 ///< parametric coordinates of the geometry

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >; ///< cell residual vector
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >; ///< cell tangent matrix

    const int _cellLabel;                                                          ///< label of the cell

    /**
     * @struct Marmot::Cells::GradientEnhancedFiniteStrainCell::MaterialPointLocation
     * @brief A hosted material point with its cached shape functions.
     */
    struct MaterialPointLocation {
      MaterialPoint*                         materialPoint; ///< the hosted material point
      XiSized                                xi;            ///< its parametric coordinates in the cell
      Fastor::Tensor< double, nNodes >       N;             ///< shape functions @f$ N_A @f$ at @c xi
      Fastor::Tensor< double, nDim, nNodes > dN_dY;         ///< gradients @f$ \partial N_A/\partial Y_i @f$
    };

    std::vector< MaterialPointLocation > _materialPointLocations; ///< the material points currently in the cell

  public:
    /**
     * @brief Construct a cell.
     * @param[in] cellLabel Label of the cell.
     * @param[in] geometry  Geometry of the cell.
     */
    GradientEnhancedFiniteStrainCell( int cellLabel, const GeometryCell& geometry )
      : GeometryCell( geometry ), _cellLabel( cellLabel ){};

    /**
     * @brief Fields per node.
     * @return For each node, `displacement` and `nonlocal damage`.
     */
    const std::vector< std::vector< std::string > >& getNodeFields() const;

    /**
     * @brief Permutation from the node-by-node dof layout of the host to the field-by-field layout of the cell.
     * @return The permutation pattern.
     */
    const std::vector< int >& getDofIndicesPermutationPattern() const;

    /**
     * @brief Supported body loads.
     * @return `BODYFORCE`.
     */
    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const { return _supportedBodyLoadTypes; }

    /**
     * @brief Supported distributed loads.
     * @return `PRESSURE`.
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
     * @brief Number of dofs of the cell.
     * @return @f$ (n_\mathrm{dim}+1)\,n_\mathrm{nodes} @f$.
     */
    int getNDofPerCell() const { return sizeLoadVector; }

    /**
     * @brief Shape of the cell.
     * @return Shape name of the geometry.
     */
    std::string getCellShape() const { return GeometryCell::getElementShape(); }

    /**
     * @brief Check whether a point lies in the cell.
     * @param[in] coordinates Coordinates of the point.
     * @return True if the point lies in the cell.
     */
    bool isCoordinateInCell( const double* coordinates ) const
    {
      return GeometryCell::isCoordinateInCell( coordinates );
    }

    /**
     * @brief Axis-aligned bounding box of the cell.
     * @param[out] boundingBoxMin Lower corner.
     * @param[out] boundingBoxMax Upper corner.
     */
    void getBoundingBox( double* boundingBoxMin, double* boundingBoxMax ) const
    {
      GeometryCell::getBoundingBox( boundingBoxMin, boundingBoxMax );
    }

    /**
     * @brief Assign the material points currently in the cell and cache their shape functions.
     *
     * The shape functions and gradients are evaluated at getCoordinatesAtCenter() of each point, i.e. at its position
     * @f$ \boldsymbol{Y} @f$ of the last accepted state.
     *
     * @param[in] materialPoints The material points.
     * @throws std::invalid_argument if a point is not a GradientEnhancedFiniteStrainMaterialPoint of this dimension.
     */
    void assignMaterialPoints( const std::vector< MarmotMaterialPoint* >& materialPoints );

    /**
     * @brief Assemble the internal force vector and its tangent from the material points.
     *
     * Requires that interpolateFieldsToMaterialPoints() and the material points' computeYourself() have been called
     * for the same increment. The contributions are ADDED to @p fInt and @p dfInt_dQ.
     *
     * @param[in]     dQ       Increment of the cell dofs (field-by-field layout).
     * @param[in,out] fInt     Internal force vector.
     * @param[in,out] dfInt_dQ Tangent, column-major.
     * @param[in]     timeNew  Time at the end of the increment (unused).
     * @param[in]     dT       Time increment (unused).
     */
    void computeMaterialPointKernels( const double* dQ,
                                      double*       fInt,
                                      double*       dfInt_dQ,
                                      double        timeNew,
                                      double        dT ) const;

    /**
     * @brief Row-sum lumped mass of the displacement field.
     * @param[out] I Lumped mass vector (overwritten); zero on the nonlocal dofs.
     */
    void computeLumpedInertia( double* I );

    /**
     * @brief Consistent mass @f$ M_{AiBi} = \sum_p N_A N_B\,\rho_0 V_p^0 @f$ of the displacement field.
     *
     * The nonlocal field carries no inertia.
     *
     * @param[out] I Mass matrix (overwritten), column-major.
     */
    void computeConsistentInertia( double* I );

    /**
     * @brief Body force load.
     *
     * `BODYFORCE`: @f$ f_{Ai} \mathrel{-}= \sum_p N_A\,b_i\,V_p^0 @f$ (sign convention of the internal force); no
     * tangent contribution.
     *
     * @param[in]     type     Load type (see getSupportedBodyLoadTypes()).
     * @param[in]     load     Body force per unit undeformed volume (nDim values).
     * @param[in,out] fExt     Load vector.
     * @param[in,out] dfExt_dQ Load tangent (untouched).
     * @param[in]     timeNew  Time at the end of the increment (unused).
     * @param[in]     dT       Time increment (unused).
     * @throws std::invalid_argument for an unknown load type.
     */
    void computeBodyLoad( int type, const double* load, double* fExt, double* dfExt_dQ, double timeNew, double dT )
      const;

    /**
     * @brief Follower pressure load at a single material point.
     *
     * `PRESSURE`: the undeformed load vector @f$ \boldsymbol{f}_0 = p\,\boldsymbol{N}\,dA_0 @f$ is pushed forward by
     * Nanson's formula, @f$ \boldsymbol{f} = J\,\boldsymbol{F}^{-\mathsf{T}}\boldsymbol{f}_0 @f$ with the total
     * @f$ \boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n @f$, and applied as
     * @f$ f_{Ai} \mathrel{-}= N_A f_i @f$, with the consistent load stiffness with respect to the displacement
     * increment.
     *
     * @param[in]     type                Load type (see getSupportedDistributedLoadTypes()).
     * @param[in]     surfaceID           Surface id (unused).
     * @param[in]     materialPointNumber Label of the loaded material point; other points are skipped.
     * @param[in]     load                Undeformed load vector @f$ \boldsymbol{f}_0 @f$ (nDim values).
     * @param[in,out] fExt                Load vector.
     * @param[in,out] dExt_dQ             Load tangent, column-major.
     * @param[in]     timeNew             Time at the end of the increment (unused).
     * @param[in]     dT                  Time increment (unused).
     * @throws std::invalid_argument for an unknown load type.
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
     * @brief Interpolate the grid increment to the material points.
     *
     * For each hosted point: @f$ \delta\boldsymbol{u} = N_A\,\Delta\boldsymbol{q}^U_A @f$,
     * @f$ \partial\delta u_i/\partial Y_j = \Delta q^U_{Ai}\,\partial N_A/\partial Y_j @f$ and
     * @f$ \delta\bar{N} = N_A\,\Delta q^N_A @f$ are passed to incrementDeformation() of the point.
     *
     * @param[in] dQ Increment of the cell dofs (field-by-field layout).
     */
    void interpolateFieldsToMaterialPoints( const double* dQ ) const;

    /**
     * @brief Cell shape functions at a point.
     * @param[out] vec         Shape functions @f$ N_A @f$ (nNodes values).
     * @param[in]  coordinates Coordinates of the point.
     */
    void getInterpolationVector( double* vec, const double* coordinates ) const;
  };

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::assignMaterialPoints(
    const std::vector< MarmotMaterialPoint* >& materialPoints )
  {
    _materialPointLocations.clear();

    XiSized coordsMP;

    for ( auto& mp : materialPoints ) {

      auto geMp = dynamic_cast< MaterialPoint* >( mp );
      if ( !geMp )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": material point " << mp->getMaterialPointNumber()
                                     << " is not a GradientEnhancedFiniteStrainMaterialPoint" );

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
  const std::vector< std::vector< std::string > >& GradientEnhancedFiniteStrainCell< nDim,
                                                                                     nNodes,
                                                                                     CellBase,
                                                                                     GeometryCell >::getNodeFields()
    const
  {
    static std::vector< std::vector< std::string > > nodeFields;

    if ( nodeFields.empty() )
      nodeFields = FiniteElement::makeNodeFieldLayout( _fields );

    return nodeFields;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  const std::vector< int >& GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::
    getDofIndicesPermutationPattern() const
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() )
      permutationPattern = FiniteElement::makeBlockedLayoutPermutationPattern( getNodeFields(), _fields );

    return permutationPattern;
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::interpolateFieldsToMaterialPoints(
    const double* dQ ) const
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;

    const auto dQU = TensorMap< const double, nNodes, nDim >( dQ );
    const auto dQN = TensorMap< const double, nNodes >( dQ + idxN );

    for ( auto& mpl : _materialPointLocations ) {

      const auto du    = evaluate( einsum< A, Ai >( mpl.N, dQU ) );
      const auto du_dY = evaluate( einsum< Ai, jA >( dQU, mpl.dN_dY ) );

      mpl.materialPoint->incrementDeformation( du, du_dY, inner( mpl.N, dQN ) );
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::computeMaterialPointKernels(
    const double* dQ,
    double*       fInt_,
    double*       dfInt_dQ_,
    double        timeNew,
    double        dT ) const
  {
    using namespace Fastor;
    using namespace Marmot::FastorIndices;

    const auto dQN = TensorMap< const double, nNodes >( dQ + idxN );

    Tensor< double, nNodes, nDim > r_U( 0.0 );
    Tensor< double, nNodes >       r_N( 0.0 );

    // the tangent is assembled directly into the output: fixed-size nNodes x nNodes blocks and their einsum temporaries
    // (~300 kB each for the displacement block of a 64-node hexahedron) overflow the 1 MB default stack of Windows
    Eigen::Map< Eigen::MatrixXd > K( dfInt_dQ_, sizeLoadVector, sizeLoadVector );

    for ( const auto& mpl : _materialPointLocations ) {

      const auto& mp = mpl.materialPoint;

      const auto& N     = mpl.N;
      const auto& dN_dY = mpl.dN_dY;

      const auto dN_dx = evaluate( einsum< ji, jA >( inv( mp->dx_dY() ), dN_dY ) );
      const auto dN_dX = evaluate( einsum< ji, jA >( mp->dY_dX(), dN_dY ) ); // independent of the increment

      const double dNonLocalField = inner( N, dQN );

      const double V0 = mp->getVolumeUndeformed();

      const auto&  S           = mp->response.S;
      const double dLocalField = mp->response.dL;
      const double c           = mp->response.nonLocalRadius * mp->response.nonLocalRadius;

      const auto& t = mp->tangents;

      // clang-format off
      const auto dS_dqU = evaluate( + einsum< ijkl, lB >( t.dS_dDeltaF, dN_dY ) );
      const auto dL_dqU = evaluate( + einsum< kl,   lB >( t.dL_dDeltaF, dN_dY ) );

      r_U  += ( + einsum< iA, ij >( dN_dx, S )                                                          ) * V0;
      r_N  += ( N * dNonLocalField + c * einsum< iA, iB, B >( dN_dX, dN_dX, dQN ) - N * dLocalField    ) * V0;
      // clang-format on

      const auto SdN_dx  = evaluate( einsum< ij, iB >( S, dN_dx ) );       // S_ij dN_B/dx_i
      const auto dSdN_dx = evaluate( einsum< ij, iA >( t.dS_dN, dN_dx ) ); // dS_ij/dN dN_A/dx_i
      const auto dNdN_dX = evaluate( einsum< iA, iB >( dN_dX, dN_dX ) );   // dN_A/dX_i dN_B/dX_i

      for ( int A = 0; A < nNodes; A++ ) {
        for ( int j = 0; j < nDim; j++ ) {
          const int rowU = idxU + A * nDim + j;
          for ( int B = 0; B < nNodes; B++ ) {
            // K_UU: ( dN_A/dx_i dS_ij/dq_Bk - dN_A/dx_k S_ij dN_B/dx_i ) V0
            for ( int k = 0; k < nDim; k++ ) {
              double kAjBk = -dN_dx( k, A ) * SdN_dx( j, B );
              for ( int i = 0; i < nDim; i++ )
                kAjBk += dN_dx( i, A ) * dS_dqU( i, j, k, B );
              K( rowU, idxU + B * nDim + k ) += kAjBk * V0;
            }
            // K_UN: dN_A/dx_i dS_ij/dN N_B V0
            K( rowU, idxN + B ) += dSdN_dx( j, A ) * N( B ) * V0;
          }
        }
        const int rowN = idxN + A;
        for ( int B = 0; B < nNodes; B++ ) {
          // K_NU: - N_A dL/dq_Bk V0
          for ( int k = 0; k < nDim; k++ )
            K( rowN, idxU + B * nDim + k ) -= N( A ) * dL_dqU( k, B ) * V0;
          // K_NN: ( N_A N_B ( 1 - dL/dN ) + c dN_A/dX_i dN_B/dX_i ) V0
          K( rowN, idxN + B ) += ( N( A ) * N( B ) * ( 1. - t.dL_dN ) + c * dNdN_dX( A, B ) ) * V0;
        }
      }
    }

    using namespace Eigen;

    // Due to Fastor Bug #139, we cannot directly write using a TensorMap
    Map< RhsSized > P( fInt_ );
    P.template segment< bsU >( idxU ) += Map< Matrix< double, bsU, 1 > >( r_U.data() );
    P.template segment< bsN >( idxN ) += Map< Matrix< double, bsN, 1 > >( r_N.data() );
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::computeConsistentInertia( double* I )
  {
    Eigen::Map< KeSizedMatrix > M( I );
    M.setZero();

    // inertia of the displacement field only; the nonlocal field is quasi-static here
    for ( const auto& mpl : _materialPointLocations ) {
      const double m = mpl.materialPoint->getDensityUndeformed() * mpl.materialPoint->getVolumeUndeformed();
      for ( int A = 0; A < nNodes; A++ )
        for ( int B = 0; B < nNodes; B++ )
          for ( int i = 0; i < nDim; i++ )
            M( idxU + A * nDim + i, idxU + B * nDim + i ) += mpl.N( A ) * mpl.N( B ) * m;
    }
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::computeLumpedInertia( double* I )
  {
    // dynamic storage: a 3D cubic B-spline cell has 256 dofs, too many for a fixed-size matrix on the stack
    Eigen::MatrixXd M( sizeLoadVector, sizeLoadVector );
    computeConsistentInertia( M.data() );
    Eigen::Map< RhsSized > lumped( I );
    lumped = M.rowwise().sum();
  }

  template < int nDim, int nNodes, class CellBase, GeometryCellPolicy< nDim, nNodes > GeometryCell >
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::computeBodyLoad( int           type,
                                                                                                  const double* load_,
                                                                                                  double*       rhs_,
                                                                                                  double* dRhs_dQ_,
                                                                                                  double  timeNew,
                                                                                                  double  dT ) const
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
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::computeDistributedLoad(
    int           type,
    int           surfaceID,
    int           materialPointNumber,
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
  void GradientEnhancedFiniteStrainCell< nDim, nNodes, CellBase, GeometryCell >::getInterpolationVector(
    double*       N,
    const double* coordinates ) const
  {
    const auto refCoord = GeometryCell::findReferenceCoordinate( XiSized( coordinates ) );

    Eigen::Map< Eigen::Matrix< double, nNodes, 1 > > interpolationVector( N );
    interpolationVector = GeometryCell::N( refCoord ).transpose();
  }

  /**
   * @class Marmot::Cells::LagrangianGradientEnhancedFiniteStrainCell
   * @brief GradientEnhancedFiniteStrainCell on a Lagrangian cell geometry (Quad4, Hexa8).
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of cell nodes.
   */
  template < int nDim, int nNodes >
  class LagrangianGradientEnhancedFiniteStrainCell
    : public GradientEnhancedFiniteStrainCell< nDim,
                                               nNodes,
                                               MarmotCell,
                                               MarmotLagrangianCellGeometry< nDim, nNodes > > {

    using Geometry      = MarmotLagrangianCellGeometry< nDim, nNodes >;                           ///< cell geometry
    using PhysicsParent = GradientEnhancedFiniteStrainCell< nDim, nNodes, MarmotCell, Geometry >; ///< physics base
                                                                                                  ///< class

  public:
    /**
     * @brief Construct a Lagrangian cell.
     * @param[in] cellLabel           Label of the cell.
     * @param[in] nodeCoordinates     Node coordinates, node by node.
     * @param[in] sizeNodeCoordinates Number of coordinates (unused).
     */
    LagrangianGradientEnhancedFiniteStrainCell( int cellLabel, const double* nodeCoordinates, int sizeNodeCoordinates )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates ) ){};
  };

  /**
   * @class Marmot::Cells::BSplineGradientEnhancedFiniteStrainCell
   * @brief GradientEnhancedFiniteStrainCell on a B-spline cell geometry.
   * @tparam nDim   Spatial dimension.
   * @tparam nNodes Number of control points of the cell, @f$ (p+1)^{n_\mathrm{dim}} @f$.
   * @tparam order  Polynomial order @f$ p @f$ of the B-splines.
   */
  template < int nDim, int nNodes, int order >
  class BSplineGradientEnhancedFiniteStrainCell
    : public GradientEnhancedFiniteStrainCell< nDim, nNodes, MarmotCell, MarmotBSplineCellGeometry< nDim, order > > {

    using Geometry      = MarmotBSplineCellGeometry< nDim, order >;                               ///< cell geometry
    using PhysicsParent = GradientEnhancedFiniteStrainCell< nDim, nNodes, MarmotCell, Geometry >; ///< physics base
                                                                                                  ///< class

  public:
    /**
     * @brief Construct a B-spline cell.
     * @param[in] cellLabel           Label of the cell.
     * @param[in] nodeCoordinates     Control point coordinates.
     * @param[in] sizeNodeCoordinates Number of coordinates.
     * @param[in] knotVectors         Knot vectors of the cell.
     * @param[in] sizeKnotVectors     Number of knot values.
     */
    BSplineGradientEnhancedFiniteStrainCell( int           cellLabel,
                                             const double* nodeCoordinates,
                                             int           sizeNodeCoordinates,
                                             const double* knotVectors,
                                             int           sizeKnotVectors )
      : PhysicsParent( cellLabel, Geometry( nodeCoordinates, sizeNodeCoordinates, knotVectors, sizeKnotVectors ) ){};
  };

} // namespace Marmot::Cells
