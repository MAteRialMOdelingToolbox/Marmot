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
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <iostream>
#include <map>

namespace Marmot::Discretizations::Elements {

  /**
   * @class GeometryElement
   * @brief Statically-sized geometry base class for all isoparametric Marmot elements.
   *
   * @tparam nDim   Number of spatial dimensions (1, 2 or 3).
   * @tparam nNodes Number of element nodes.
   *
   * Provides shape functions @c N, their natural-coordinate derivatives @c dNdXi,
   * the B-operator @c B, and the element Jacobian.  The element shape is determined
   * automatically from the template parameters.
   *
   * Concrete Marmot::Discretizations::Elements::Element classes inherit from this class to access the geometric
   * infrastructure without code duplication.
   */
  template < int nDim, int nNodes >
  class GeometryElement {

  public:
    /*Typedefs*/
    /* static constexpr VoigtSize voigtSize = ( ( ( nDim * nDim ) + nDim ) / 2 ); */
    /// @brief Voigt notation size for @p nDim spatial dimensions.
    static constexpr Marmot::ContinuumMechanics::Voigt::VoigtSize
      voigtSize = Marmot::ContinuumMechanics::Voigt::voigtSizeFromDimension( nDim );

    typedef Eigen::Matrix< double, nDim, 1 >             XiSized;          ///< Natural-coordinate vector.
    typedef Eigen::Matrix< double, nDim * nNodes, 1 >    CoordinateVector; ///< Flat nodal-coordinate vector.
    typedef Eigen::Matrix< double, nDim, nDim >          JacobianSized;    ///< Square Jacobian matrix.
    typedef Eigen::Matrix< double, 1, nNodes >           NSized;           ///< Row vector of shape function values.
    typedef Eigen::Matrix< double, nDim, nNodes * nDim > NBSized;          ///< Expanded interpolation operator N_B.
    typedef Eigen::Matrix< double, nDim, nNodes >        dNdXiSized;  ///< Matrix of shape function natural derivatives.
    typedef Eigen::Matrix< double, voigtSize, nNodes * nDim > BSized; ///< Standard B-operator matrix.
    typedef Eigen::Matrix< double, 4, nNodes * nDim >
      BSizedAxisymmetric; ///< Axisymmetric B-operator matrix (4 Voigt components).

    /*Properties*/
    Eigen::Map< const CoordinateVector > coordinates; ///< Map into the externally-owned nodal-coordinate array.

    /// @brief Element shape determined from @p nDim and @p nNodes, at compile time.
    static constexpr ElementShapes shape = getElementShapeByMetric( nDim, nNodes );
    /// @brief Reference domain of the shape: the cube [-1, 1]^nDim, the line [-1, 1] or the unit simplex.
    static constexpr ReferenceDomain referenceDomain = getReferenceDomain( shape );
    /// @brief Polynomial order of the shape functions.
    static constexpr int interpolationOrder = getInterpolationOrder( shape );

    /*Methods*/
    /// @brief Default constructor. Initialises the coordinate map to @c nullptr.
    GeometryElement() : coordinates( nullptr ){};

    /// @brief Returns an Ensight Gold shape string (e.g. @c quad4, @c hexa8) for this element.
    std::string getElementShape() const
    {
      using namespace Marmot::Discretizations::Elements;
      static std::map< ElementShapes, std::string > shapes = { { Bar2, "bar2" },
                                                               { Tri3, "tria3" },
                                                               { Quad4, "quad4" },
                                                               { Quad8, "quad8" },
                                                               { Tetra4, "tetra4" },
                                                               { Tetra10, "tetra10" },
                                                               { Hexa8, "hexa8" },
                                                               { Hexa20, "hexa20" } };

      return shapes[this->shape];
    }

    /// @brief Maps the externally-owned coordinate array into the @ref coordinates member.
    /// @param[in] coords Pointer to the flat nodal-coordinate array.
    void assignNodeCoordinates( const double* coords )
    {
      new ( &coordinates ) Eigen::Map< const CoordinateVector >( coords );
    }

    /*Please specialize these functions for each element individially
     *.cpp file.
     *Fully specialized templates are precompiled in marmotMechanics (rather than the unspecialized and
     *partially specialized templates)
     * */
    /// @brief Evaluate the shape function row vector at natural coordinates @p xi.
    NSized N( const XiSized& xi ) const;
    /// @brief Evaluate the matrix of shape function natural derivatives at @p xi.
    dNdXiSized dNdXi( const XiSized& xi ) const;
    /// @brief Compute the standard B-operator from physical derivatives @p dNdX.
    BSized B( const dNdXiSized& dNdX ) const;
    /// @brief Compute the axisymmetric B-operator (4 components).
    BSizedAxisymmetric B_axisymmetric( const dNdXiSized& dNdX, const NSized& N, const XiSized& x_gauss ) const;
    /// @brief Compute the B̄-operator using the B-bar selective-reduced-integration method.
    BSized B_bar( const dNdXiSized& dNdX, const dNdXiSized& dNdX0 ) const;
    /// @brief Compute the Green–Lagrange strain operator for deformation gradient @p F.
    BSized BGreen( const dNdXiSized& dNdX, const JacobianSized& F ) const;

    /*These functions are equal for each element and independent of node number and  nDimension*/
    /// @brief Compute the expanded interpolation operator N_B from shape function values @p N.
    NBSized NB( const NSized& N ) const { return Marmot::Discretizations::Elements::NB< nDim, nNodes >( N ); }

    /// @brief Compute the element Jacobian from natural derivatives @p dNdXi and the stored nodal coordinates.
    JacobianSized Jacobian( const dNdXiSized& dNdXi ) const
    {
      return Marmot::Discretizations::Elements::Jacobian< nDim, nNodes >( dNdXi, coordinates );
    }

    /// @brief Compute physical derivatives @p dN/dX from natural derivatives and the inverse Jacobian.
    /// @param[in] dNdXi        Natural-coordinate derivatives.
    /// @param[in] JacobianInverse  Inverse of the element Jacobian.
    dNdXiSized dNdX( const dNdXiSized& dNdXi, const JacobianSized& JacobianInverse ) const
    {
      return ( dNdXi.transpose() * JacobianInverse ).transpose();
    }

    /// @brief Compute the deformation gradient \f$\mathbf{F}\f$ from physical derivatives and displacements @p Q.
    /// @param[in] dNdX Physical shape function derivatives.
    /// @param[in] Q    Element displacement vector.
    JacobianSized F( const dNdXiSized& dNdX, const CoordinateVector& Q ) const
    {
      return Marmot::Discretizations::Elements::Jacobian< nDim, nNodes >( dNdX, Q ) + JacobianSized::Identity();
    }

    /*Element geometry depending on the reference domain: the only place that distinguishes the cube from the simplex,
     * so the elements do not have to*/

    /// @brief The natural coordinates of the centroid of the reference domain: the origin of the cube [-1, 1]^nDim,
    /// ( 1 / ( nDim + 1 ), ... ) for the unit simplex.
    static XiSized referenceCentroid()
    {
      return XiSized::Constant( getReferenceCentroidCoordinate( referenceDomain ) );
    }

    /// @brief The physical coordinates of the element center, the image of the centroid of the reference domain.
    /// @details For a linear simplex all shape functions are 1 / nNodes there, so it is the mean of the nodal
    /// coordinates. Evaluating at the origin instead would give node 0 of a simplex. It is the true centroid (the
    /// center of area or volume) only for an affine element, i.e. a linear simplex, a parallelogram or a
    /// parallelepiped; for a general quadrilateral or hexahedron it is the mean of the corner nodes.
    Eigen::Matrix< double, nDim, 1 > center() const { return NB( N( referenceCentroid() ) ) * coordinates; }

    /**
     * @brief The characteristic element length \f$L_e\f$ at the natural coordinates @p xi, for the stable time
     * increment \f$\Delta t = f\, L_e / c\f$ of explicit dynamics, with the wave speed \f$c\f$ and a factor \f$f\f$
     * for the lumped mass distribution (see MassLumping::timeStepFactorFromMassDistribution).
     *
     * @details
     * - Quadrilateral, hexahedron and line: \f$L_e = 2\,\sigma_\mathrm{min}(\mathbf{J})\f$, twice the smallest
     *   singular value of the Jacobian. The Jacobian maps the cube \f$[-1,1]^{n_\mathrm{dim}}\f$ onto the element, so
     *   this is the element's smallest physical extent.
     * - Linear simplex: \f$L_e = 1 / \sqrt{\sum_i |\nabla N_i|^2} = 1 / \sqrt{\sum_i h_i^{-2}}\f$, with the constant
     *   shape function gradients \f$\nabla N_i\f$ and the altitude \f$h_i\f$ opposite node \f$i\f$. Unlike the
     *   Jacobian of the unit simplex, it does not depend on the node numbering.
     *
     * Quadratic simplices (Tri6, Tetra10) have no derived length; using one is a compile error.
     *
     * @param[in] xi Natural coordinates; irrelevant for a linear simplex, whose gradients are constant.
     */
    double characteristicLengthForCriticalTimeStepAt( const XiSized& xi ) const
    {
      static_assert( !isSimplex( referenceDomain ) || interpolationOrder == 1,
                     "The characteristic element length of a quadratic simplex is not derived." );

      const dNdXiSized    dNdXi_ = dNdXi( xi );
      const JacobianSized J      = Jacobian( dNdXi_ );

      if constexpr ( isSimplex( referenceDomain ) ) {
        // the squared Frobenius norm of the nDim x nNodes gradient matrix is the sum over i of |grad N_i|^2
        return 1.0 / std::sqrt( dNdX( dNdXi_, J.inverse() ).squaredNorm() );
      }
      else {
        // A Jacobian with non-finite entries, from invalid nodal coordinates, leaves the singular values unset
        // (InvalidInput), so they are read only after a successful decomposition.
        const Eigen::JacobiSVD< JacobianSized > svd( J );
        if ( svd.info() != Eigen::Success )
          throw std::invalid_argument(
            "GeometryElement::characteristicLengthForCriticalTimeStepAt(): the Jacobian is not finite; "
            "check the nodal coordinates." );
        return 2.0 * svd.singularValues().minCoeff();
      }
    }

    /**
     * @brief The volume-equivalent length \f$(V_\mathrm{ref} \det\mathbf{J})^{1/n_\mathrm{dim}}\f$: the side of the
     * square or cube with the measure the element would have if its Jacobian were @p detJ everywhere.
     *
     * @details \f$V_\mathrm{ref}\f$ is the measure of the reference domain (getReferenceMeasure()): 2, 4 and 8 for the
     * line and the cubes \f$[-1,1]^{n_\mathrm{dim}}\f$, 1/2 and 1/6 for the unit triangle and tetrahedron, for
     * elements of any order on that domain. For an element with straight edges \f$\det\mathbf{J}\f$ is constant and
     * this is exactly \f$\sqrt{A}\f$ or \f$\sqrt[3]{V}\f$; for a curved or distorted one it is the measure the element
     * would have with the Jacobian at that point. Elements pass it to their materials as the regularisation length of
     * softening laws.
     *
     * @param[in] detJ Determinant of the Jacobian, e.g. at a quadrature point.
     */
    static double volumeEquivalentLength( double detJ )
    {
      const double measure = getReferenceMeasure( referenceDomain ) * detJ;
      if constexpr ( nDim == 3 )
        return std::cbrt( measure );
      else if constexpr ( nDim == 2 )
        return std::sqrt( measure );
      else
        return measure;
    }
  };

} // namespace Marmot::Discretizations::Elements
