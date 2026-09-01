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
 * Alexandros Stathas alexandros.stathas@boku.ac.at
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
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

/**
 * @class MarmotGeometryInterfaceElement
 * @brief Geometry helper for zero-thickness interface elements.
 *
 * @tparam nDim Spatial embedding dimension of the interface element. Supported values are 2 and 3.
 * @tparam nNodes Total number of element nodes. The first half are the bottom side, the second half are the top side.
 *
 * @details
 * The class evaluates shape functions, surface Jacobians, metric tensors, normals, surface gradients, and the
 * interpolation operators required by interface finite elements. For an interface embedded in `nDim` dimensions, the
 * parametric interface dimension is `nDim - 1`.
 *
 * Supported interpolations are:
 * - `nDim == 2`, `nNodes == 4`: two-node line interface (`iline2`).
 * - `nDim == 3`, `nNodes == 8`: four-node quadrilateral interface (`iquad4`).
 *
 * Node ordering is expected to be side-wise:
 * - nodes `0 ... nInterfaceNodes - 1`: bottom side,
 * - nodes `nInterfaceNodes ... nNodes - 1`: top side.
 */
template < int nDim, int nNodes >
class MarmotGeometryInterfaceElement {
public:
  static_assert( nDim == 2 || nDim == 3,
                 "MarmotGeometryInterfaceElement supports interfaces embedded in 2D or 3D only." );

  static_assert( nNodes % 2 == 0, "Interface element must have bottom and top sides, so nNodes must be even." );

  /// @brief Number of parametric coordinates on the interface manifold.
  static constexpr int nXi = nDim - 1;

  /// @brief Number of interpolation nodes per interface side.
  static constexpr int nInterfaceNodes = nNodes / 2;

  /// @brief Number of displacement degrees of freedom per node.
  static constexpr int nDofPerNodeU = nDim;

  /// @brief Number of displacement degrees of freedom on one interface side.
  static constexpr int nDofPerSide = nInterfaceNodes * nDim;

  /// @brief Number of displacement degrees of freedom of the full interface element.
  static constexpr int nDofElement = nNodes * nDim;

  /// @brief Parametric coordinate vector on the interface.
  using XiSized = Eigen::Matrix< double, nXi, 1 >;

  /// @brief Row vector of scalar shape functions for one interface side.
  using NSized = Eigen::Matrix< double, 1, nInterfaceNodes >;

  /// @brief Derivatives of scalar shape functions with respect to parametric coordinates.
  using dNdXiSized = Eigen::Matrix< double, nXi, nInterfaceNodes >;

  /// @brief Stacked coordinate vector of all interface element nodes.
  using CoordinateVector = Eigen::Matrix< double, nDim * nNodes, 1 >;

  /// @brief Stacked coordinate vector of one interface side.
  using SideCoordinateVector = Eigen::Matrix< double, nDim * nInterfaceNodes, 1 >;

  /// @brief Surface Jacobian mapping parametric interface directions into physical space.
  using SurfaceJacobianSized = Eigen::Matrix< double, nDim, nXi >;

  /// @brief First fundamental form, `G = J.transpose() * J`.
  using MetricSized = Eigen::Matrix< double, nXi, nXi >;

  /// @brief Surface gradients of scalar shape functions in physical coordinates.
  using GradSized = Eigen::Matrix< double, nDim, nInterfaceNodes >;

  /// @brief Vector in the embedding space.
  using VectorDim = Eigen::Matrix< double, nDim, 1 >;

  /// @brief Second-order tensor in the embedding space.
  using TensorDim = Eigen::Matrix< double, nDim, nDim >;

  /// @brief Vector-valued interpolation matrix for one interface side.
  using NMatrixSized = Eigen::Matrix< double, nDim, nDofPerSide >;

  /// @brief Interpolation matrix for the displacement jump `u_top - u_bottom`.
  using NJumpMatrixSized = Eigen::Matrix< double, nDim, nDofElement >;

  /// @brief Surface displacement-gradient matrix for one interface side.
  using BSurfaceSized = Eigen::Matrix< double, nDim * nDim, nDofPerSide >;

  /// @brief Surface displacement-gradient matrix for the average of bottom and top sides.
  using BAvgSurfaceSized = Eigen::Matrix< double, nDim * nDim, nDofElement >;

  /// @brief Non-owning map onto the stacked element coordinates.
  Eigen::Map< const CoordinateVector > coordinates;

  /// @brief Internal Marmot shape identifier for the computational interface interpolation.
  const Marmot::FiniteElement::ElementShapes shape;

  /**
   * @brief Construct an interface geometry object without assigned coordinates.
   *
   * @details
   * Coordinates must be attached later with assignNodeCoordinates() before any geometry evaluation is performed.
   */
  MarmotGeometryInterfaceElement() : coordinates( nullptr ), shape( getInterfaceElementShape() ) {}

  /**
   * @brief Return the Marmot element shape associated with the template parameters.
   *
   * @return `Bar2` for 2D line interfaces and `Quad4` for 3D quadrilateral interfaces.
   */
  static constexpr Marmot::FiniteElement::ElementShapes getInterfaceElementShape()
  {
    using namespace Marmot::FiniteElement;

    if constexpr ( nDim == 2 && nInterfaceNodes == 2 ) {
      return Bar2;
    }
    else if constexpr ( nDim == 3 && nInterfaceNodes == 4 ) {
      return Quad4;
    }
    else {
      static_assert( nDim == -1, "Unsupported MarmotGeometryInterfaceElement interpolation in the current build." );
    }
  }

  /**
   * @brief Return the result-file geometry keyword for this interface interpolation.
   *
   * @return `"iline2"` for 2D line interfaces and `"iquad4"` for 3D quadrilateral interfaces.
   */
  std::string getElementShape() const
  {
    using namespace Marmot::FiniteElement;

    static std::map< ElementShapes, std::string > shapes = {
      { Bar2, "iline2" },
      { Quad4, "iquad4" },
    };

    return shapes[this->shape];
  }

  /**
   * @brief Attach nodal coordinates to this geometry object.
   *
   * @param coords Pointer to `nDim * nNodes` contiguous coordinate values in side-wise node order.
   *
   * @warning The coordinate memory is not copied and must remain valid for the lifetime of this geometry object or
   * until assignNodeCoordinates() is called again.
   */
  void assignNodeCoordinates( const double* coords )
  {
    new ( &coordinates ) Eigen::Map< const CoordinateVector >( coords );
  }

  /**
   * @brief Evaluate scalar shape functions at a parametric interface point.
   *
   * @param xi Parametric coordinates on the interface.
   * @return Row vector containing one shape function value per interface-side node.
   *
   * @note Specializations are provided for supported interface interpolations.
   */
  NSized N( const XiSized& xi ) const;

  /**
   * @brief Evaluate parametric derivatives of scalar shape functions.
   *
   * @param xi Parametric coordinates on the interface.
   * @return Matrix of derivatives `dN_A / dxi_alpha`.
   *
   * @note Specializations are provided for supported interface interpolations.
   */
  dNdXiSized dNdXi( const XiSized& xi ) const;

  /**
   * @brief Extract the physical coordinates of one interface side.
   *
   * @param side Interface side index, where `0` is the bottom side and `1` is the top side.
   * @return Stacked coordinate vector of the selected side.
   *
   * @throws std::invalid_argument If `side` is neither `0` nor `1`.
   */
  SideCoordinateVector getSideCoordinates( const int side = 0 ) const
  {
    if ( side != 0 && side != 1 )
      throw std::invalid_argument( "Interface side must be 0 or 1." );

    SideCoordinateVector xSide;

    const int nodeOffset = side * nInterfaceNodes;

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        xSide( A * nDim + i ) = coordinates( ( nodeOffset + A ) * nDim + i );
      }
    }

    return xSide;
  }

  /**
   * @brief Compute the surface Jacobian of one interface side.
   *
   * @param dN Parametric shape-function derivatives.
   * @param side Interface side index used for the geometry, where `0` is bottom and `1` is top.
   * @return Surface Jacobian `J`, with columns spanning the physical tangent space.
   */
  SurfaceJacobianSized surfaceJacobian( const dNdXiSized& dN, const int side = 0 ) const
  {
    const SideCoordinateVector xSide = getSideCoordinates( side );

    SurfaceJacobianSized J;
    J.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        const double XAi = xSide( A * nDim + i );

        for ( int alpha = 0; alpha < nXi; ++alpha ) {
          J( i, alpha ) += XAi * dN( alpha, A );
        }
      }
    }

    return J;
  }

  /**
   * @brief Compute the metric tensor associated with a surface Jacobian.
   *
   * @param J Surface Jacobian.
   * @return Metric tensor `G = J.transpose() * J`.
   */
  MetricSized metric( const SurfaceJacobianSized& J ) const { return J.transpose() * J; }

  /**
   * @brief Compute the square root of the metric determinant.
   *
   * @param J Surface Jacobian.
   * @return Surface measure factor `sqrt(det(G))`.
   *
   * @throws std::runtime_error If the metric determinant is non-positive.
   */
  double sqrtDetMetric( const SurfaceJacobianSized& J ) const
  {
    const MetricSized G    = metric( J );
    const double      detG = G.determinant();

    if ( detG <= 0.0 )
      throw std::runtime_error( "Degenerate interface element: det(G) <= 0." );

    return std::sqrt( detG );
  }

  /**
   * @brief Compute physical surface gradients of scalar shape functions.
   *
   * @param dN Parametric shape-function derivatives.
   * @param J Surface Jacobian.
   * @return Matrix whose columns are `grad_s N_A`.
   */
  GradSized surfaceGradient( const dNdXiSized& dN, const SurfaceJacobianSized& J ) const
  {
    const MetricSized GInv = metric( J ).inverse();

    /*
     * Surface gradient of scalar shape functions:
     *
     *   grad_s N_A = J * G^{-1} * dN_A/dxi
     *
     * Shape:
     *   J       : nDim x nXi
     *   GInv    : nXi  x nXi
     *   dN      : nXi  x nInterfaceNodes
     *   gradN   : nDim x nInterfaceNodes
     */
    return J * GInv * dN;
  }

  /**
   * @brief Compute the unit normal vector of the interface side.
   *
   * @param J Surface Jacobian.
   * @return Unit normal vector. In 2D this is a 90-degree rotation of the tangent; in 3D it is the normalized cross
   * product of the two tangent vectors.
   *
   * @throws std::runtime_error If the normal vector has near-zero norm.
   */
  VectorDim normal( const SurfaceJacobianSized& J ) const
  {
    VectorDim n;
    n.setZero();

    if constexpr ( nDim == 2 ) {
      const VectorDim t = J.col( 0 );

      n( 0 ) = -t( 1 );
      n( 1 ) = t( 0 );
    }
    else if constexpr ( nDim == 3 ) {
      n = J.col( 0 ).cross( J.col( 1 ) );
    }

    const double normN = n.norm();

    if ( normN < 1e-16 )
      throw std::runtime_error( "Degenerate interface element: normal vector is zero." );

    return n / normN;
  }

  /**
   * @brief Compute the normal projection tensor.
   *
   * @param n Unit normal vector.
   * @return Projection tensor `n * n.transpose()`.
   */
  TensorDim normalProjector( const VectorDim& n ) const { return n * n.transpose(); }

  /**
   * @brief Compute the tangent projection tensor.
   *
   * @param n Unit normal vector.
   * @return Projection tensor `I - n * n.transpose()`.
   */
  TensorDim tangentProjector( const VectorDim& n ) const { return TensorDim::Identity() - normalProjector( n ); }

  /**
   * @brief Build the vector-valued interpolation matrix for one interface side.
   *
   * @param N_ Scalar shape-function row vector.
   * @return Matrix mapping side nodal displacements to interpolated displacement.
   */
  NMatrixSized NMatrix( const NSized& N_ ) const
  {
    NMatrixSized Nmat;
    Nmat.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        Nmat( i, A * nDim + i ) = N_( 0, A );
      }
    }

    return Nmat;
  }

  /**
   * @brief Build the displacement-jump interpolation matrix.
   *
   * @param N_ Scalar shape-function row vector.
   * @return Matrix mapping full element nodal displacements to `u_top - u_bottom`.
   */
  NJumpMatrixSized NJumpMatrix( const NSized& N_ ) const
  {
    NJumpMatrixSized NJump;
    NJump.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int i = 0; i < nDim; ++i ) {
        const int colBottom = A * nDim + i;
        const int colTop    = ( nInterfaceNodes + A ) * nDim + i;

        NJump( i, colBottom ) = -N_( 0, A );
        NJump( i, colTop )    = N_( 0, A );
      }
    }

    return NJump;
  }

  /**
   * @brief Build the projected surface displacement-gradient matrix for one side.
   *
   * @param gradN Physical surface gradients of scalar shape functions.
   * @param T Tangent projection tensor.
   * @return Matrix mapping side nodal displacements to the projected full surface displacement gradient.
   *
   * @details
   * The operator stores the surface displacement gradient in row-major tensor order
   * `row = i * nDim + k`, where `i` is the displacement component and `k` is the gradient direction. It is not the
   * symmetric small-strain `B` matrix.
   *
   * Only the GRADIENT DIRECTION is projected: `B(i,k),(A,i) = gradN(j,A) T(j,k)`. The displacement component `i`
   * is left unprojected, so the operator produces the surface gradient of the FULL displacement vector,
   * `grad_s u = grad(u) . T`, and retains the in-plane derivative of the normal component.
   *
   * That is the measure MarmotInterfaceMaterialHypoElastic expects, since it forms the thin-layer gradient as
   * `grad(u) = (1/h) [u] (x) n + <grad_s u>`, in which the normal displacement varying along the surface is a
   * genuine layer shear. Projecting the displacement component as well (`T(i,m) gradN(j,A) T(j,k)`) would delete
   * that term. The same formula is used for every `nDim`, so a plane-strain model and its 3D extrusion produce
   * the same surface strain.
   */
  BSurfaceSized BSurfaceMatrix( const GradSized& gradN, const TensorDim& T ) const
  {
    BSurfaceSized B;
    B.setZero();

    for ( int A = 0; A < nInterfaceNodes; ++A ) {
      for ( int k = 0; k < nDim; ++k ) {

        double value = 0.0;

        for ( int j = 0; j < nDim; ++j ) {
          value += gradN( j, A ) * T( j, k );
        }

        for ( int i = 0; i < nDim; ++i ) {
          const int row = i * nDim + k;
          const int col = A * nDim + i;

          B( row, col ) = value;
        }
      }
    }

    return B;
  }

  /**
   * @brief Build the average surface displacement-gradient matrix for the full interface element.
   *
   * @param Bside Surface displacement-gradient matrix for one interface side.
   * @return Matrix applying `0.5 * Bside` to both bottom and top side degrees of freedom.
   */
  BAvgSurfaceSized BAverageSurfaceMatrix( const BSurfaceSized& Bside ) const
  {
    BAvgSurfaceSized BAvg;
    BAvg.setZero();

    for ( int row = 0; row < nDim * nDim; ++row ) {
      for ( int col = 0; col < nDofPerSide; ++col ) {
        BAvg( row, col )               = 0.5 * Bside( row, col );
        BAvg( row, nDofPerSide + col ) = 0.5 * Bside( row, col );
      }
    }

    return BAvg;
  }

  /**
   * @struct QuadratureGeometry
   * @brief Bundle of geometry quantities evaluated at one interface quadrature point.
   *
   * @details
   * The struct stores scalar shape functions, differential geometry terms, projection tensors, and interpolation
   * operators needed by interface element assembly. All members are value types with compile-time dimensions.
   */
  struct QuadratureGeometry {
    /// @brief Scalar shape-function row vector.
    NSized N;

    /// @brief Parametric derivatives of scalar shape functions.
    dNdXiSized dNdXi;

    /// @brief Surface Jacobian of the side selected for geometry evaluation.
    SurfaceJacobianSized J;

    /// @brief Metric tensor associated with `J`.
    MetricSized G;

    /// @brief Surface measure factor `sqrt(det(G))`.
    double sqrtDetG;

    /// @brief Physical surface gradients of scalar shape functions.
    GradSized gradN;

    /// @brief Unit normal vector.
    VectorDim n;

    /// @brief Normal projection tensor `n * n.transpose()`.
    TensorDim normalProjection;

    /// @brief Tangent projection tensor `I - n * n.transpose()`.
    TensorDim tangentProjection;

    /// @brief Vector-valued interpolation matrix for one interface side.
    NMatrixSized NmatSide;

    /// @brief Displacement-jump interpolation matrix for the full interface element.
    NJumpMatrixSized NmatJump;

    /// @brief Surface displacement-gradient matrix for one interface side.
    BSurfaceSized BmatSide;

    /// @brief Average surface displacement-gradient matrix for the full interface element.
    BAvgSurfaceSized BmatAverage;
  };

  /**
   * @brief Evaluate all geometry quantities needed at one quadrature point.
   *
   * @param xi Parametric coordinates of the quadrature point.
   * @param sideForGeometry Interface side used for the surface Jacobian and derived geometry terms. Use `0` for
   * bottom-side geometry and `1` for top-side geometry.
   * @return Fully populated QuadratureGeometry bundle.
   *
   * @throws std::invalid_argument If `sideForGeometry` is neither `0` nor `1`.
   * @throws std::runtime_error If the selected side is geometrically degenerate.
   */
  QuadratureGeometry evaluateAt( const XiSized& xi, const int sideForGeometry = 0 ) const
  {
    QuadratureGeometry q;

    q.N     = N( xi );
    q.dNdXi = dNdXi( xi );

    q.J        = surfaceJacobian( q.dNdXi, sideForGeometry );
    q.G        = metric( q.J );
    q.sqrtDetG = sqrtDetMetric( q.J );

    q.gradN = surfaceGradient( q.dNdXi, q.J );

    q.n = normal( q.J );

    q.normalProjection  = normalProjector( q.n );
    q.tangentProjection = tangentProjector( q.n );

    q.NmatSide = NMatrix( q.N );
    q.NmatJump = NJumpMatrix( q.N );

    q.BmatSide = BSurfaceMatrix( q.gradN, q.tangentProjection );

    q.BmatAverage = BAverageSurfaceMatrix( q.BmatSide );

    return q;
  }
};
