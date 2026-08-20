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
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include <Fastor/tensor_algebra/permute.h>

namespace Marmot::ContinuumMechanics {

  namespace DeformationMeasures {

    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    /** @brief Computes the right Cauchy-Green tensor \f$\boldsymbol{C}\f$.
     * @tparam T Scalar type (e.g., float, double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return The right Cauchy-Green tensor
     *
     * The right Cauchy-Green deformation tensor is defined as:
     * \f[
     *   \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F}
     * \f]
     * or in index notation:
     * \f[
     *   C_{IJ} = F_{iI} F_{iJ}
     * \f]
     */
    template < typename T >
    Tensor33t< T > rightCauchyGreen( const Tensor33t< T >& F )
    {
      const Tensor33t< T > C = einsum< iI, iJ >( F, F );
      return C;
    }

    /** @brief Computes the left Cauchy-Green tensor \f$\boldsymbol{b}\f$.
     * @tparam T Scalar type (e.g., float, double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return The left Cauchy-Green tensor b
     *
     * The left Cauchy-Green deformation tensor is defined as:
     * \f[
     *   \boldsymbol{b} = \boldsymbol{F} \boldsymbol{F}^T
     * \f]
     * or in index notation:
     * \f[
     *   b_{ij} = F_{iI} F_{jI}
     * \f]
     */
    template < typename T >
    Tensor33t< T > leftCauchyGreen( const Tensor33t< T >& F )
    {
      const Tensor33t< T > b = einsum< iJ, jJ >( F, F );
      return b;
    }

    /** @brief Computes the principal stretches \f$\lambda_i\f$ of \f$\boldsymbol{F}\f$ and the
     * corresponding eigenvector matrix \f$\boldsymbol{Q}\f$ of the right Cauchy-Green tensor
     * \f$\boldsymbol{C}\f$ (the material principal frame).
     * @tparam T Scalar type (e.g., double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return A pair containing the principal stretches \f$\lambda_i = \sqrt{\mathrm{eig}_i(C)}\f$ and
     * the eigenvector matrix \f$\boldsymbol{Q}\f$ of \f$\boldsymbol{C}\f$.
     *
     * This factors out the spectral-decomposition step shared by every principal-stretch-based
     * hyperelastic material (e.g. isotropic Biot-type materials, Ogden-type materials): decompose
     * \f$\boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F}\f$ into eigenvalues and eigenvectors via
     * Marmot::Math::computeEigenSystemJacobi(), then take the square root of the eigenvalues to
     * obtain the principal stretches. Callers assemble whatever principal-basis quantity they need
     * (e.g. the right stretch tensor \f$\boldsymbol{U} = \boldsymbol{Q}\,\mathrm{diag}(\lambda_i)\,
     * \boldsymbol{Q}^T\f$) from the returned stretches and eigenvectors.
     */
    template < typename T >
    std::pair< Tensor3t< T >, Tensor33t< T > > principalStretches( const Tensor33t< T >& F )
    {
      const Tensor33t< T > C = rightCauchyGreen( F );
      auto [eigenvalues, Q]  = Marmot::Math::computeEigenSystemJacobi( C );
      Tensor3t< T > principalStretch;
      for ( int i = 0; i < 3; ++i )
        principalStretch( i ) = sqrt( eigenvalues( i ) );
      return { principalStretch, Q };
    }

    /** @brief Computes the volume ratio (Jacobian determinant) \f$J = \det(\boldsymbol{F})\f$.
     * @tparam T Scalar type (e.g., double, autodiff::dual)
     * @param F Deformation gradient tensor
     * @return The volume ratio J
     */
    template < typename T >
    T volumeRatio( const Tensor33t< T >& F )
    {
      return determinant( F );
    }

    /** @brief Computes the F-bar modified deformation gradient
     * \f$\bar{\boldsymbol{F}} = (J_0/J)^{1/3}\,\boldsymbol{F}\f$ (de Souza Neto et al., 1996).
     *
     * @details Replacing \f$\boldsymbol{F}\f$ by \f$\bar{\boldsymbol{F}}\f$ before evaluating a
     * material's constitutive response decouples the volumetric response (sampled once via
     * \f$J_0\f$, e.g. at the element centroid) from the deviatoric response (sampled at each
     * quadrature point via \f$J\f$), which alleviates volumetric locking of nearly incompressible
     * materials on low-order displacement elements (e.g. Marmot::Elements::
     * DisplacementFiniteStrainFBarElement). This is the value only; the element-level consistent
     * tangent additionally needs a correction term coupling \f$J_0\f$'s dependence on every
     * node's degrees of freedom, which is not a per-quadrature-point quantity and is therefore
     * assembled in the element itself rather than here.
     *
     * @tparam T Scalar type (e.g., double, autodiff::dual)
     * @param F Deformation gradient tensor at the current quadrature point.
     * @param J Volume ratio \f$J=\det\boldsymbol F\f$ at the current quadrature point.
     * @param J0 Volume ratio sampled at the element centroid (or another single reference point).
     * @return The F-bar modified deformation gradient \f$\bar{\boldsymbol{F}}\f$.
     */
    template < typename T >
    Tensor33t< T > FBarDeformationGradient( const Tensor33t< T >& F, const T& J, const T& J0 )
    {
      return pow( J0 / J, 1. / 3. ) * F;
    }

    namespace FirstOrderDerived {

      /** @brief Computes the right Cauchy-Green tensor \f$\boldsymbol{C}\f$ and its derivative with respect to
       * \f$\boldsymbol{F}\f$.
       * @tparam T Scalar type (e.g., float, double, autodiff::dual)
       * @param F Deformation gradient tensor
       * @returns A pair containing the right Cauchy-Green tensor C and its derivative with respect to F
       *
       * The right Cauchy-Green deformation tensor is defined as:
       * \f[
       *   \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F}
       * \f]
       * or in index notation:
       * \f[
       *   C_{IJ} = F_{iI} F_{iJ}
       * \f]
       * The derivative of \f$\boldsymbol{C}\f$ with respect to \f$\boldsymbol{F}\f$ is given by:
       * \f[
       *  \frac{\partial C_{IJ}}{\partial F_{kK}} = \delta_{IK} F_{kJ} + F_{kI} \delta_{JK}
       * \f]
       */
      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > rightCauchyGreen( const Tensor33t< T >& F )
      {
        const Tensor33t< T > C = Marmot::ContinuumMechanics::DeformationMeasures::rightCauchyGreen( F );

        const Tensor3333t< T > dC_dF = einsum< IK, kJ, to_IJkK >( Spatial3D::I, F ) +
                                       einsum< kI, JK, to_IJkK >( F, Spatial3D::I );

        return { C, dC_dF };
      }

      /** @brief Computes the left Cauchy-Green tensor \f$\boldsymbol{b}\f$ and its derivative with respect to
       * \f$\boldsymbol{F}\f$.
       * @tparam T Scalar type (e.g., float, double, autodiff::dual)
       * @param F Deformation gradient tensor
       * @return A pair containing the left Cauchy-Green tensor b and its derivative with respect to F
       *
       * The left Cauchy-Green deformation tensor is defined as:
       * \f[
       *   \boldsymbol{b} = \boldsymbol{F} \boldsymbol{F}^T
       * \f]
       * or in index notation:
       * \f[
       *   b_{ij} = F_{iI} F_{jI}
       * \f]
       * The derivative of \f$\boldsymbol{b}\f$ with respect to \f$\boldsymbol{F}\f$ is given by:
       * \f[
       *  \frac{\partial b_{ij}}{\partial F_{kK}} = \delta_{ik} F_{jK} + F_{iK} \delta_{jk}
       * \f]
       */
      template < typename T >
      std::pair< Tensor33t< T >, Tensor3333t< T > > leftCauchyGreen( const Tensor33t< T >& F )
      {
        Tensor33t< T > b = Marmot::ContinuumMechanics::DeformationMeasures::leftCauchyGreen( F );

        const Tensor3333t< T > db_dF = einsum< ik, jK, to_ijkK >( Spatial3D::I, F ) +
                                       einsum< iK, jk, to_ijkK >( F, Spatial3D::I );

        return { b, db_dF };
      }

      /** @brief Computes the deformed normal vector and its derivative with respect to the inverse
       * deformation gradient.
       * @param FInv Inverse of the deformation gradient tensor
       * @param N_x_dA0 Undeformed normal vector scaled by undeformed area
       * @return A pair containing the deformed normal vector and its derivative with respect to FInv
       *
       * The deformed normal vector is defined as:
       * \f[
       *   n_i = \frac{ F^{-1}_{Ii} N_{I} J dA_0 } { \| F^{-1}_{Ki} N_{K} J dA_0 \| }
       * \f]
       */

      template < int nDim >
      std::pair< Fastor::Tensor< double, nDim >, Fastor::Tensor< double, nDim, nDim, nDim > > deformedNormalVectorFromUndeformedSurfaceVector(
        const Fastor::Tensor< double, nDim, nDim >& FInv,
        const Fastor::Tensor< double, nDim >&       N_x_dA0 )
      {

        using namespace Fastor;
        using namespace FastorIndices;

        using TensorDd   = Tensor< double, nDim >;
        using TensorDDd  = Tensor< double, nDim, nDim >;
        using TensorDDDd = Tensor< double, nDim, nDim, nDim >;

        const static TensorDDd Eye = []() {
          TensorDDd I;
          I.eye();
          return I;
        }();

        const double    detJ        = 1. / determinant( FInv );
        const TensorDDd dDetJ_dFInv = -1 * std::pow( detJ, 2 ) *
                                      transpose( inverse( ( FInv ) ) ); // derivative of detJ w.r.t. FInv

        const TensorDd
                         NTilde        = einsum< Index< 0 >, Index< 0, 1 > >( N_x_dA0,
                                                        FInv ); // deformed normal vector scaled by deformed area
        const TensorDDDd dNTilde_dFInv = einsum< Index< 0 >, Index< 1, 2 >, OIndex< 1, 0, 2 > >( N_x_dA0, Eye );

        const TensorDd   n_x_dA               = detJ * NTilde;
        const TensorDDDd NTilde_x_dDetJ_dFInv = outer( NTilde, dDetJ_dFInv );
        const TensorDDDd dn_x_dA_dFInv        = NTilde_x_dDetJ_dFInv + detJ * dNTilde_dFInv;

        const TensorDd  n          = n_x_dA / norm( n_x_dA ); // deformed normal vector
        const TensorDDd dn_dn_x_dA = ( Eye - outer( n, n ) ) * ( 1. / norm( n_x_dA ) );

        const TensorDDDd dn_dFInv = einsum< ij, jkl >( dn_dn_x_dA, dn_x_dA_dFInv );

        return { n, dn_dFInv };
      }

      /** @brief Computes the deformed normal projection tensor and its derivative with respect to the inverse
       * deformation gradient.
       * @param FInv Inverse of the deformation gradient tensor
       * @param N_x_dA0 Undeformed normal vector scaled by undeformed area
       * @return A pair containing the deformed normal projection tensor and its derivative with respect to FInv
       *
       * The deformed normal projection tensor is defined as:
       * \f[
       *   N_{ij} = n_i n_j = \frac{ F^{-1}_{Ii} N_{I} J dA_0 \times F^{-1}_{Jj} N_{J} J dA_0} { \| F^{-1}_{Ki} N_{K} J
       * dA_0 \|^2 }
       * \f]
       */
      template < int nDim >
      inline std::pair< Fastor::Tensor< double, nDim, nDim >, Fastor::Tensor< double, nDim, nDim, nDim, nDim > > deformedNormalProjectionTensor(
        const Fastor::Tensor< double, nDim, nDim >& FInv,
        const Fastor::Tensor< double, nDim >&       N_x_dA0 )
      {

        using namespace Fastor;
        using namespace FastorIndices;
        using TensorDDd   = Tensor< double, nDim, nDim >;
        using TensorDDDDd = Tensor< double, nDim, nDim, nDim, nDim >;

        const auto [n, dn_dFInv] = deformedNormalVectorFromUndeformedSurfaceVector< nDim >( FInv, N_x_dA0 );
        const TensorDDd n_ij     = outer( n, n );

        const TensorDDDDd aux         = outer( n, dn_dFInv );
        const TensorDDDDd dn_ij_dFInv = aux + Fastor::permute< Fastor::Index< 1, 0, 2, 3 > >( aux );

        return { n_ij, dn_ij_dFInv };
      }

      /**
       * @brief Computes the inverse of the deformation gradient and its derivative with respect to the deformation
       * gradient.
       * @param F Deformation gradient tensor
       * @return A pair containing the inverse of the deformation gradient and its derivative with respect to F
       *
       * The inverse of the deformation gradient is defined as:
       * \f[
       *  \boldsymbol{F}^{-1}
       *  \f]
       *  The derivative of \f$\boldsymbol{F}^{-1}\f$ with respect to \f$\boldsymbol{F}\f$ is given by:
       *  \f[
       *  \frac{\partial F^{-1}_{iK}}{\partial F_{jL}} = - F^{-1}_{iL} F^{-1}_{jK}
       *  \f]
       */
      template < int nDim >
      inline std::pair< Fastor::Tensor< double, nDim, nDim >, Fastor::Tensor< double, nDim, nDim, nDim, nDim > > inverseDeformationGradient(
        const Fastor::Tensor< double, nDim, nDim >& F )
      {

        using namespace Fastor;
        using namespace FastorIndices;
        using TensorDDd   = Tensor< double, nDim, nDim >;
        using TensorDDDDd = Tensor< double, nDim, nDim, nDim, nDim >;

        const TensorDDd   FInv     = inverse( F );
        const TensorDDDDd dFInv_dF = -einsum< Ik, Ki, to_IikK >( FInv, FInv );

        return { FInv, dFInv_dF };
      }

      /** @brief Computes the volume ratio \f$J = \det(\boldsymbol{F})\f$ and its derivative with
       * respect to \f$\boldsymbol{F}\f$.
       * @tparam T Scalar type (e.g., double, autodiff::dual)
       * @param F Deformation gradient tensor
       * @return A pair containing \f$J\f$ and \f$\partial J/\partial\boldsymbol{F} = J\,
       * \boldsymbol{F}^{-T}\f$.
       */
      template < typename T >
      std::pair< T, Tensor33t< T > > volumeRatio( const Tensor33t< T >& F )
      {
        const T              J     = determinant( F );
        const Tensor33t< T > FInv  = inverse( F );
        const Tensor33t< T > dJ_dF = multiplyFastorTensorWithScalar( transpose( FInv ), J );

        return { J, dJ_dF };
      }

    } // namespace FirstOrderDerived

  }   // namespace DeformationMeasures
} // namespace Marmot::ContinuumMechanics
