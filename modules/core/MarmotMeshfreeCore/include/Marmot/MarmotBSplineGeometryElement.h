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
#include "Eigen/Core"
#include "Marmot/MarmotBSpline.h"
#include "Marmot/MarmotFiniteElement.h"
#include <iostream>
#include <map>

namespace Marmot::FiniteElement {

  /**
   * @class Marmot::FiniteElement::MarmotBSplineGeometryElement
   * @brief Geometry of a single knot span of a tensor-product B-spline of degree @p order.
   *
   * @details The element is one knot span @f$ [u_p, u_{p+1}] @f$ per direction of a B-spline of degree
   * @f$ p @f$ = @p order. The @f$ p+1 @f$ basis functions that are nonzero on that span depend on the
   * @f$ 2p+2 @f$ knots @f$ u_0,\dots,u_{2p+1} @f$ of each direction, which are passed to the constructor. The
   * shape functions are tensor products of the univariate B-splines @f$ B_{i,p} @f$ (Cox--de Boor recursion,
   * see MarmotBSpline.h),
   * @f[
   *   N_{a}(\boldsymbol{\xi}) = B_{i,p}(\xi_1)\,B_{j,p}(\xi_2)\,B_{k,p}(\xi_3), \qquad
   *   a = i + j\,(p+1) + k\,(p+1)^2 ,
   * @f]
   * i.e. with the index of the first direction running fastest. The parametric coordinates @f$ \boldsymbol{\xi} @f$
   * are the knot values themselves (the normalization of the knot vectors to @f$ [-1,1] @f$ is disabled, so
   * #_knotVectorsNormalized is a copy of #_knotVectors).
   *
   * @note Only 2D and 3D are usable: in 1D, the constructor throws for @p order > 1 (no shape with @f$ p+1 > 2 @f$
   * nodes in 1D), and dNdXi() would not compile if instantiated.
   *
   * @tparam nDim  Spatial dimension (2 or 3).
   * @tparam order Polynomial degree @f$ p @f$ of the B-spline.
   */
  template < int nDim, int order >
  class MarmotBSplineGeometryElement {

    /**
     * @brief Compile-time integer power.
     * @param[in] base     Base.
     * @param[in] exponent Non-negative exponent.
     * @return @f$ \mathrm{base}^{\mathrm{exponent}} @f$.
     */
    static int constexpr pow( int base, int exponent ) { return exponent == 0 ? 1 : base * pow( base, exponent - 1 ); }

  public:
    constexpr static int nKnotsPerDir = 2 * order + 2; ///< Knots per direction, @f$ 2p+2 @f$.
    using KnotVector                  = Eigen::Matrix< double, nKnotsPerDir, 1 >; ///< Knot vector of one direction.
    using KnotVectors = Eigen::Matrix< double, nKnotsPerDir, nDim >; ///< Knot vectors, one column per direction.

    constexpr static int nNodes = pow( order + 1, nDim ); ///< Number of control points, @f$ (p+1)^{n_\mathrm{dim}} @f$.
    using NSized                = Eigen::Matrix< double, 1, nNodes >;        ///< Row vector of shape function values.
    using CoordinateVector      = Eigen::Matrix< double, nDim * nNodes, 1 >; ///< Flat control point coordinates.
    using JacobianSized         = Eigen::Matrix< double, nDim, nDim >;       ///< Jacobian matrix.

    using dNdXSized = Eigen::Matrix< double, nDim, nNodes >;    ///< Shape function derivatives (nDim x nNodes).
    using XiSized   = Eigen::Matrix< double, nDim, 1 >;         ///< Parametric point.

    const Eigen::Map< const CoordinateVector > _mapCoordinates; ///< Map to the externally owned control point
                                                                ///< coordinates (point by point).
    const Marmot::FiniteElement::ElementShapes _shape;          ///< Shape deduced from nDim and nNodes.

    KnotVectors _knotVectors;                                   ///< Knot vectors as passed to the constructor.
    KnotVectors _knotVectorsNormalized; ///< Knot vectors used for evaluation; currently a copy of #_knotVectors.

    /**
     * @brief Constructs the geometry of one knot span.
     * @param[in] nodeCoordinates_    Control point coordinates, nNodes @f$ \times @f$ nDim values, point by point in
     *                                the tensor-product order of N(). Mapped, not copied; must outlive the object.
     * @param[in] sizeNodeCoordinates Size of @p nodeCoordinates_ (not checked).
     * @param[in] knotVectors_        Knot vectors, nKnotsPerDir values per direction, direction by direction.
     * @param[in] sizeKnotVectors     Size of @p knotVectors_ (not checked).
     * @throws std::invalid_argument if no element shape exists for (nDim, nNodes).
     */
    MarmotBSplineGeometryElement( const double* nodeCoordinates_,
                                  int           sizeNodeCoordinates,
                                  const double* knotVectors_,
                                  int           sizeKnotVectors )
      : _mapCoordinates( Eigen::Map< const CoordinateVector >( nodeCoordinates_ ) ),
        _shape( Marmot::FiniteElement::getElementShapeByMetric( nDim, nNodes ) )
    {
      _knotVectors = Eigen::Map< const KnotVectors >( knotVectors_ );

      for ( int d = 0; d < nDim; d++ ) {

        /* knotVectors_ */
        _knotVectorsNormalized = _knotVectors;
        /* _knotVectorsNormalized.col( */
        /*   d ) = ( _knotVectors.array().col( d ) - _knotVectors.col( d )( order ) ) * */
        /*         ( 2 / ( _knotVectors.col( d )( nKnotsPerDir - order - 1 ) - _knotVectors.col( d )( order ) ) ); */

        /* _knotVectorsNormalized.array().col( d ) -= 1; */
      }
    };

    /**
     * @brief Ensight Gold shape name used for output.
     * @return @c "quad4", @c "quad9", @c "hexa8" or @c "hexa27" for degree 1 and 2; an empty string for degree 3
     *         (@c Quad16 / @c Hexa64 have no entry).
     */
    std::string getElementShape() const
    {
      using namespace Marmot::FiniteElement;
      static std::map< ElementShapes, std::string > shapes = {
        { Bar2, "bar2" },
        { Quad4, "quad4" },
        { Quad8, "quad8" },
        { Quad9, "quad9" },
        { Hexa8, "hexa8" },
        { Hexa20, "hexa20" },
        { Hexa27, "hexa27" },
      };

      return shapes[this->_shape];
    }

    /**
     * @brief Tensor-product B-spline shape functions.
     * @param[in] xi Parametric coordinates (knot values), inside the knot span.
     * @return @f$ N_a(\boldsymbol{\xi}) @f$.
     */
    NSized N( const XiSized& xi ) const
    {
      NSized        N_;
      constexpr int nN = order + 1;

      if constexpr ( nDim == 1 )
        for ( int p = 0; p < nN; p++ )
          N_( p ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p );

      else if constexpr ( nDim == 2 )
        for ( int q = 0; q < nN; q++ )
          for ( int p = 0; p < nN; p++ )
            N_( p + q * ( nN ) ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                   B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );

      else if constexpr ( nDim == 3 )
        for ( int r = 0; r < nN; r++ )
          for ( int q = 0; q < nN; q++ )
            for ( int p = 0; p < nN; p++ )
              N_( p + q * nN + r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                               B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                               B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );

      return N_;
    }

    /**
     * @brief Derivatives of the shape functions w.r.t. the parametric coordinates.
     * @param[in] xi Parametric coordinates (knot values), inside the knot span.
     * @return @f$ \partial N_a/\partial \xi_i @f$ (row @f$ i @f$, column @f$ a @f$).
     */
    dNdXSized dNdXi( const XiSized& xi ) const
    {
      dNdXSized     dN_dXi_;
      constexpr int nN = order + 1;

      if constexpr ( nDim == 1 )
        for ( int p = 0; p < nN; p++ )
          dN_dXi_( p ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized[0].data(), p );

      else if constexpr ( nDim == 2 )
        for ( int q = 0; q < nN; q++ )
          for ( int p = 0; p < nN; p++ ) {
            dN_dXi_( 0, p + q * nN ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                       B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );
            dN_dXi_( 1, p + q * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                       dB_dU< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q );
          }

      else if constexpr ( nDim == 3 )
        for ( int r = 0; r < nN; r++ )
          for ( int q = 0; q < nN; q++ )
            for ( int p = 0; p < nN; p++ ) {
              dN_dXi_( 0,
                       p + q * nN +
                         r * nN * nN ) = dB_dU< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
              dN_dXi_( 1,
                       p + q * nN +
                         r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         dB_dU< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         B< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
              dN_dXi_( 2,
                       p + q * nN +
                         r * nN * nN ) = B< order >( xi( 0 ), this->_knotVectorsNormalized.col( 0 ).data(), p ) *
                                         B< order >( xi( 1 ), this->_knotVectorsNormalized.col( 1 ).data(), q ) *
                                         dB_dU< order >( xi( 2 ), this->_knotVectorsNormalized.col( 2 ).data(), r );
            }

      return dN_dXi_;
    }

    /**
     * @brief Jacobian of the geometry map defined by the control points.
     * @param[in] dNdXi Shape function derivatives from dNdXi().
     * @return @f$ J_{ij} = \sum_a \partial N_a/\partial \xi_j \, X_{a,i} @f$.
     */
    JacobianSized Jacobian( const dNdXSized& dNdXi ) const
    {
      return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdXi, _mapCoordinates );
    }

    /**
     * @brief Transforms parametric to physical shape function derivatives.
     * @param[in] dNdXi           Shape function derivatives from dNdXi().
     * @param[in] JacobianInverse Inverse of Jacobian().
     * @return @f$ \partial N_a/\partial X_i = (\partial N_a/\partial \xi_j)\,J^{-1}_{ji} @f$.
     */
    dNdXSized dNdX( const dNdXSized& dNdXi, const JacobianSized& JacobianInverse ) const
    {
      return ( dNdXi.transpose() * JacobianInverse ).transpose();
    }
  };

} // namespace Marmot::FiniteElement
