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
#include <iostream>
#include <map>

namespace Marmot {

  /**
   * @class MarmotGeometryElement
   * @brief Statically-sized geometry base class for all isoparametric Marmot elements.
   *
   * @tparam nDim   Number of spatial dimensions (1, 2 or 3).
   * @tparam nNodes Number of element nodes.
   *
   * Provides shape functions @c N, their natural-coordinate derivatives @c dNdXi,
   * the B-operator @c B, and the element Jacobian.  The element shape is determined
   * automatically from the template parameters.
   *
   * Concrete MarmotElement classes inherit from this class to access the geometric
   * infrastructure without code duplication.
   */
  template < int nDim, int nNodes >
  class MarmotGeometryElement {

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
    Eigen::Map< const CoordinateVector >       coordinates; ///< Map into the externally-owned nodal-coordinate array.
    const Marmot::FiniteElement::ElementShapes shape;       ///< Element shape determined from @p nDim and @p nNodes.

    /*Methods*/
    /// @brief Default constructor. Initialises the coordinate map to @c nullptr and deduces the element shape.
    MarmotGeometryElement()
      : coordinates( nullptr ), shape( Marmot::FiniteElement::getElementShapeByMetric( nDim, nNodes ) ){};

    /// @brief Returns an Ensight Gold shape string (e.g. @c quad4, @c hexa8) for this element.
    std::string getElementShape() const
    {
      using namespace FiniteElement;
      static std::map< ElementShapes, std::string > shapes = { { Bar2, "bar2" },
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
    NBSized NB( const NSized& N ) const { return Marmot::FiniteElement::NB< nDim, nNodes >( N ); }

    /// @brief Compute the element Jacobian from natural derivatives @p dNdXi and the stored nodal coordinates.
    JacobianSized Jacobian( const dNdXiSized& dNdXi ) const
    {
      return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdXi, coordinates );
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
      return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdX, Q ) + JacobianSized::Identity();
    }
  };

} // namespace Marmot
