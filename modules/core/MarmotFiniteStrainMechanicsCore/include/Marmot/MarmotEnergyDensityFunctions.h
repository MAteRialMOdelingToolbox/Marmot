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
#include "Marmot/MarmotFastorTensorBasics.h"

namespace Marmot::ContinuumMechanics {

  namespace EnergyDensityFunctions {

    using namespace Fastor;
    using namespace FastorStandardTensors;

    /**
     * @brief Hyperelastic potentials that include an explicit volumetric penalty \f$U(J)\f$ and
     * are intended to be driven by a plain displacement element.
     */
    namespace Compressible {

      /** @brief Hyperelastic Energy Density Function Wa acc. Pence & Gou (2015), Eq. (2.11)
       *
       *  The energy density function \f$W_a\f$ is given as
       *  \f[
       *    W_a = \frac{G}{2} (I_1 - 3) + \left(\frac{K}{2} - \frac{G}{3}\right) (J - 1)^2 - G \ln(J)
       *  \f]
       *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
       *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
       * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
       * shear modulus, respectively.
       *
       * @tparam T Scalar type, e.g. double, float, etc.
       * @param C Right Cauchy-Green tensor
       * @param K Bulk modulus
       * @param G Shear modulus
       * @return Energy density
       */
      template < typename T >
      T PenceGouPotentialA( const Tensor33t< T >& C, const double K, const double G )
      {

        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );

        T res = G / 2. * ( I1 - 3. ) + ( K / 2. - G / 3. ) * pow( J - 1, 2 ) - G * log( J );

        return res;
      }

      /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12)
       *
       *  The energy density function \f$W_b\f$ is given as
       *  \f[
       *    W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
       *  \f]
       *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
       *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
       * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
       * shear modulus, respectively.
       *
       * @tparam T Scalar type, e.g. double, float, etc.
       * @param C Right Cauchy-Green tensor
       * @param K Bulk modulus
       * @param G Shear modulus
       * @return Energy density
       */
      template < typename T >
      T PenceGouPotentialB( const Tensor33t< T >& C, const double K, const double G )
      {

        const T detC = determinant( C );
        const T I1   = trace( C );

        T res = K / 8. * ( detC + 1. / detC - 2. ) + G / 2. * ( I1 * pow( detC, -1. / 3 ) - 3. );

        return res;
      }

      /** @brief Hyperelastic Energy Density Function Wc acc. Pence & Gou (2015), Eq. (2.13)
       *
       *  The energy density function \f$W_c\f$ is given as
       *  \f[
       *    W_c = \frac{G}{2} (I_1 - 3) + \frac{3 G^2}{3 K - 2 G} \left(J^{\frac{2}{3} - \frac{K}{G}} - 1\right)
       *  \f]
       *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
       *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
       * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
       * shear modulus, respectively.
       *
       * @tparam T Scalar type, e.g. double, float, etc.
       * @param C Right Cauchy-Green tensor
       * @param K Bulk modulus
       * @param G Shear modulus
       * @return Energy density
       */
      template < typename T >
      T PenceGouPotentialC( const Tensor33t< T >& C, const double K, const double G )
      {

        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );

        T res = G / 2. * ( I1 - 3. ) + 3. * G * G / ( 3. * K - 2. * G ) * ( pow( J, 2. / 3 - K / G ) - 1 );

        return res;
      }

      /**
       * @brief Mooney-Rivlin Hyperelastic Energy Density Function
       * The energy density function \f$W_{MR}\f$ is given as
       * \f[
       *   W_{MR} = C_1 (\bar{I}_1 - 3) + C_2 (\bar{I}_2 - 3) + \frac{1}{D_1} (J - 1)^2
       * \f]
       * where \f$ \bar{I}_1 = I_1 J^{-\frac{2}{3}} \f$ and \f$ \bar{I}_2 = I_2 J^{-\frac{4}{3}} \f$ are the first
       * and second invariant of the isochoric right Cauchy-Green tensor, respectively, \f$ I_1 =
       * \text{tr}(\boldsymbol{C}) \f$ and \f$ I_2 = 0.5 (I_1^2 - \text{tr}(\boldsymbol{C}^2)) \f$ are the first and
       * second invariant of the right Cauchy-Green tensor \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$
       * J = \sqrt{\det(\boldsymbol{C})} = \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and
       * \f$ C_1, C_2, D_1 \f$ are material parameters.
       *
       * @tparam T Scalar type, e.g. double, float, etc.
       * @param C Right Cauchy-Green tensor
       * @param C1 Mooney-Rivlin material parameter C1
       * @param C2 Mooney-Rivlin material parameter C2
       * @param D1 Mooney-Rivlin material parameter D1
       * @return Energy density
       */
      template < typename T >
      T MooneyRivlinPotential( const Tensor33t< T >& C, const double C1, const double C2, const double D1 )
      {

        const T J   = sqrt( determinant( C ) );
        const T I1  = trace( C );
        const T I1_ = I1 * pow( J, -2. / 3. );
        const T I2_ = 0.5 * ( I1 * I1 - trace( C % C ) ) * pow( J, -4. / 3. );
        T       res = C1 * ( I1_ - 3. ) + C2 * ( I2_ - 3. ) + 1. / D1 * ( 0.5 * ( J * J - 1 ) - log( J ) );

        return res;
      }

      /** @brief Yeoh hyperelastic energy density function.
       *
       *  The energy density function is given as
       *  \f[
       *    W = C_1 (\bar{I}_1 - 3) + C_2 (\bar{I}_1 - 3)^2 + C_3 (\bar{I}_1 - 3)^3 + \frac{1}{D_1}\left(
       *    \frac{J^2 - 1}{2} - \ln J \right)
       *  \f]
       *  where \f$\bar{I}_1 = I_1 J^{-2/3}\f$, \f$I_1 = \mathrm{tr}(\boldsymbol{C})\f$, and
       *  \f$J = \sqrt{\det(\boldsymbol{C})}\f$.
       *
       * @tparam T Scalar type, e.g. double, float, autodiff scalar.
       * @param C Right Cauchy-Green tensor.
       * @param C1 Yeoh material parameter.
       * @param C2 Yeoh material parameter.
       * @param C3 Yeoh material parameter.
       * @param D1 Volumetric penalty parameter.
       * @return Energy density.
       */
      template < typename T >
      T YeohPotential( const Tensor33t< T >& C, const double C1, const double C2, const double C3, const double D1 )
      {

        const T J         = sqrt( determinant( C ) );
        const T I1        = trace( C );
        const T I1_minus3 = I1 * pow( J, -2. / 3. ) - 3.;
        T       res       = C1 * I1_minus3 + C2 * I1_minus3 * I1_minus3 + C3 * I1_minus3 * I1_minus3 * I1_minus3 +
                1. / D1 * ( 0.5 * ( J * J - 1 ) - log( J ) );
        return res;
      }

      /** @brief Standard compressible Neo-Hooke energy density function in terms of \f$\boldsymbol{C}\f$.
       *
       * @tparam T Scalar type, e.g. double, float, autodiff scalar.
       * @param C Right Cauchy-Green tensor.
       * @param K Bulk modulus.
       * @param G Shear modulus.
       * @return Energy density.
       */
      template < typename T >
      T standardNeoHooke( const Tensor33t< T >& C, const double K, const double G )
      {
        const double lambda = K - 2.0 / 3.0 * G;

        /*
         * we use the potential in terms of C
         * Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + lambda/2 ( 0.5 ( J^2 - 1 ) - ln(J) )
         * where J = sqrt( det(C) ) and ln(J) = 0.5 ln( det(C) )
         *
         * we can the rewrite as
         * Psi = G/2 ( tr(C) - 3 - ln(det(C)) ) + lambda/2 ( 0.5 ( det(C) - 1 ) - 0.5 ln(det(C)) )
         *
         */

        const T trC    = trace( C );
        const T detC   = determinant( C );
        const T lnDetC = log( detC );
        // energy density
        const T psi = G / 2 * ( trC - 3.0 - lnDetC ) + lambda / 4 * ( detC - 1.0 - lnDetC );

        return psi;
      }

      namespace FirstOrderDerived {

        /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12) and its first derivative
         * w.r.t. C
         *
         *  The energy density function \f$W_b\f$ is given as
         *  \f[
         *    W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
         *  \f]
         *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
         *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
         * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
         * shear modulus, respectively.
         *
         *  Additionally, the first derivative w.r.t. C is computed as
         *  \f[
         *    \frac{\partial W_b}{\partial \boldsymbol{C}} = \frac{\partial W_b}{\partial J} \frac{\partial J}{\partial
         * \boldsymbol{C}} + \frac{\partial W_b}{\partial I_1} \frac{\partial I_1}{\partial \boldsymbol{C}} \f] where
         * \f[
         *    \frac{\partial J}{\partial \boldsymbol{C}} = \frac{1}{2} J \boldsymbol{C}^{-1}
         *  \f]
         *  and
         *  \f[
         *    \frac{\partial I_1}{\partial \boldsymbol{C}} = \boldsymbol{I}
         *  \f]
         *
         * @tparam T Scalar type, e.g. double, float, etc.
         * @param C Right Cauchy-Green tensor
         * @param K Bulk modulus
         * @param G Shear modulus
         * @return A tuple containing energy density and its first derivative w.r.t. C
         */
        template < typename T >
        std::tuple< T, Tensor33t< T > > PenceGouPotentialB( const Tensor33t< T >& C, const double K, const double G )
        {
          using namespace FastorIndices;

          const T J  = sqrt( determinant( C ) );
          const T I1 = trace( C );
          // energy density
          T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

          // first derivative w.r.t. C
          const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1. + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
          const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

          const Tensor33t< T > CInv   = inverse( C );
          const Tensor33t< T > dJ_dC  = multiplyFastorTensorWithScalar( transpose( CInv ), T( J / 2. ) );
          const Tensor33t< T > dI1_dC = fastorTensorFromDoubleTensor< T >( Spatial3D::I );

          Tensor33t< T > dPsi_dC = multiplyFastorTensorWithScalar( dJ_dC, dPsi_dJ ) +
                                   multiplyFastorTensorWithScalar( dI1_dC, dPsi_dI1 );

          return { psi, dPsi_dC };
        }
      } // namespace FirstOrderDerived

      namespace SecondOrderDerived {

        /** @brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12) and its first and second
         * derivative w.r.t. C
         *
         * The energy density function \f$W_b\f$ is given as
         * \f[
         *   W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
         * \f]
         * where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
         * \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
         * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
         * shear modulus, respectively.
         *
         * Additionally, the first and second derivative w.r.t. \f$\boldsymbol{C}\f$,
         * i.e.,
         * \f[
         * \frac{\partial W_b}{\partial \boldsymbol{C}} \quad \text{and} \quad
         * \frac{\partial^2 W_b}{\partial \boldsymbol{C} \partial \boldsymbol{C}}
         * \f]
         * are computed.
         *
         * @tparam T Scalar type, e.g. double, float, etc.
         * @param C Right Cauchy-Green tensor
         * @param K Bulk modulus
         * @param G Shear modulus
         * @return A tuple containing energy density, its first and second derivative w.r.t. C
         *
         */
        template < typename T >
        std::tuple< T, Tensor33t< T >, Tensor3333t< T > > PenceGouPotentialB( const Tensor33t< T >& C,
                                                                              const double          K,
                                                                              const double          G )
        {
          using namespace FastorIndices;

          const T J  = sqrt( determinant( C ) );
          const T I1 = trace( C );
          // energy density
          T psi = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

          // first derivative w.r.t. C
          const T dPsi_dJ  = K / 4. * ( J - 1. / J ) * ( 1 + 1. / ( J * J ) ) - G / 3. * I1 * pow( J, -5. / 3. );
          const T dPsi_dI1 = G / 2. * pow( J, -2. / 3. );

          const Tensor33t< T > CInv   = inverse( C );
          const Tensor33t< T > dJ_dC  = 0.5 * J * transpose( CInv );
          const Tensor33t< T > dI1_dC = Spatial3D::I;

          Tensor33t< T > dPsi_dC = dPsi_dJ * dJ_dC + dPsi_dI1 * dI1_dC;

          // second derivative w.r.t. C
          const T d2Psi_dJdJ        = K / 4. * ( 1. + 3. / ( J * J * J * J ) ) + 5. / 9. * G * I1 * pow( J, -8. / 3. );
          const T d2Psi_dJdI1       = -G / 3. * pow( J, -5. / 3. );
          Tensor3333t< T > d2J_dCdC = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                      J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );

          Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) + dPsi_dJ * d2J_dCdC +
                                        d2Psi_dJdI1 *
                                          ( einsum< IJ, KL >( dJ_dC, dI1_dC ) + einsum< IJ, KL >( dI1_dC, dJ_dC ) );

          return { psi, dPsi_dC, d2Psi_dCdC };
        }

        /** @brief Standard compressible Neo-Hooke energy density and first/second derivatives w.r.t.
         * \f$\boldsymbol{C}\f$.
         *
         * @tparam T Scalar type, e.g. double, float, autodiff scalar.
         * @param C Right Cauchy-Green tensor.
         * @param K Bulk modulus.
         * @param G Shear modulus.
         * @return Tuple of energy density, first derivative and second derivative w.r.t. \f$\boldsymbol{C}\f$.
         */
        template < typename T >
        std::tuple< T, FastorStandardTensors::Tensor33t< T >, FastorStandardTensors::Tensor3333t< T > > standardNeoHooke(
          const FastorStandardTensors::Tensor33t< T >& C,
          const double&                                K,
          const double&                                G )
        {
          const double lambda = K - 2.0 / 3.0 * G;

          /*
           * we use the potential in terms of C
           * Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + lambda/2 ( 0.5 ( J^2 - 1 ) - ln(J) )
           * where J = sqrt( det(C) ) and ln(J) = 0.5 ln( det(C) )
           *
           * we can the rewrite as
           * Psi = G/2 ( tr(C) - 3 - ln(det(C)) ) + lambda/2 ( 0.5 ( det(C) - 1 ) - 0.5 ln(det(C)) )
           *
           */

          const T trC    = trace( C );
          const T detC   = determinant( C );
          const T lnDetC = log( detC );
          // energy density
          const T psi = G / 2 * ( trC - 3.0 - lnDetC ) + lambda / 4 * ( detC - 1.0 - lnDetC );

          // first derivative quantities
          const Tensor33t< T >& dTrC_dC    = makeOtherScalarType< T >( Spatial3D::I );
          const Tensor33t< T >  invCt      = transpose( inverse( C ) );
          const Tensor33t< T >  dDetC_dC   = multiplyFastorTensorWithScalar( invCt, detC );
          const Tensor33t< T >& dLnDetC_dC = invCt;

          // first derivative with respect to C
          const Tensor33t< T > dPsi_dC = G / 2 * ( dTrC_dC - dLnDetC_dC ) + lambda / 4 * ( dDetC_dC - dLnDetC_dC );

          // second derivative quantities
          using namespace FastorIndices;
          using lj                          = Index< l_, j_ >;
          using li                          = Index< l_, i_ >;
          const Tensor3333t< T > dInvCt_dC  = -0.5 * ( einsum< ik, lj, to_ijkl >( invCt, invCt ) +
                                                      einsum< jk, li, to_ijkl >( invCt, invCt ) );
          const Tensor3333t< T > d2DetC_dC2 = multiplyFastorTensorWithScalar( dInvCt_dC, detC ) +
                                              einsum< ij, kl, to_ijkl >( invCt, dDetC_dC );
          const Tensor3333t< T >& d2LnDetC_dC2 = dInvCt_dC;

          // second derivative with respect to C
          const Tensor3333t< T > d2Psi_dC2 = lambda / 4.0 * d2DetC_dC2 - ( G / 2.0 + lambda / 4.0 ) * d2LnDetC_dC2;
          return { psi, dPsi_dC, d2Psi_dC2 };
        }

        /** @brief Isotropic Biot-Neo-Hooke energy density and derivatives w.r.t. right stretch \f$\boldsymbol{U}\f$.
         *
         * @details Internally evaluates the standard Neo-Hooke potential in terms of
         * \f$\boldsymbol{C}=\boldsymbol{U}\boldsymbol{U}\f$ and applies chain-rule transformations.
         *
         * @tparam T Scalar type, e.g. double, float, autodiff scalar.
         * @param U Right stretch tensor.
         * @param K Bulk modulus.
         * @param G Shear modulus.
         * @return Tuple of energy density, first derivative and second derivative w.r.t. \f$\boldsymbol{U}\f$.
         */
        template < typename T >
        std::tuple< T, FastorStandardTensors::Tensor33t< T >, FastorStandardTensors::Tensor3333t< T > > BiotNeoHooke(
          const FastorStandardTensors::Tensor33t< T >& U,
          const double&                                K,
          const double&                                G )
        {
          using namespace FastorIndices;
          Tensor3333t< T > I4       = makeOtherScalarType< T >( Spatial3D::I4 );
          Tensor33t< T >   C        = U % U;
          Tensor3333t< T > dC_dU    = 2 * einsum< iL, ijkl >( U, I4 );
          Tensor3333t< T > d2C_dUdU = 2 * I4;

          auto [psi, dPsi_dC, d2Psi_dCdC] = standardNeoHooke< T >( C, K, G );

          Tensor33t< T >   dPsi_dU    = einsum< kl, klmn >( dPsi_dC, dC_dU );
          Tensor3333t< T > d2Psi_dUdU = einsum< ijkl, klmn >( einsum< ijkl, ijmn >( dC_dU, d2Psi_dCdC ), dC_dU ) +
                                        einsum< jL, ijkl >( dPsi_dC, d2C_dUdU );

          return { psi, dPsi_dU, d2Psi_dUdU };
        }

      } // namespace SecondOrderDerived

      namespace ThirdOrderDerived {
        /** @brief Standard compressible Neo-Hooke energy density and first/second/third derivatives w.r.t.
         * \f$\boldsymbol{C}\f$.
         *
         * @param C Right Cauchy-Green tensor.
         * @param K Bulk modulus.
         * @param G Shear modulus.
         * @return Tuple of energy density, first derivative, second derivative and third derivative.
         */
        std::tuple< double, FastorStandardTensors::Tensor33d, FastorStandardTensors::Tensor3333d, FastorStandardTensors::Tensor333333d > standardNeoHooke(
          const FastorStandardTensors::Tensor33d& C,
          const double&                           K,
          const double&                           G );
      } // namespace ThirdOrderDerived

    }   // namespace Compressible

    /**
     * @brief Purely isochoric hyperelastic potentials with no volumetric term, meant to be paired
     * with a mixed displacement/pressure element that enforces incompressibility as an independent
     * constraint (see e.g. Marmot::Elements::DisplacementPressureFiniteStrainElement).
     */
    namespace Incompressible {

      /** @brief Classical N-term Ogden isochoric hyperelastic potential and its principal
       * Kirchhoff stress.
       *
       * The isochoric energy density is
       * \f[
       *   \Psi_\mathrm{iso}(\bar\lambda) = \sum_{p=1}^N \frac{\mu_p}{\alpha_p} \left(
       *   \bar\lambda_1^{\alpha_p} + \bar\lambda_2^{\alpha_p} + \bar\lambda_3^{\alpha_p} - 3 \right)
       * \f]
       * with isochoric principal stretches \f$\bar\lambda_i = J^{-1/3}\lambda_i\f$, and \f$J =
       * \det(\boldsymbol{F})\f$. There is no volumetric term \f$U(J)\f$: this potential is
       * purely isochoric and must be paired with an independent pressure field (a mixed
       * displacement-pressure element) to enforce incompressibility.
       *
       * The principal Kirchhoff stress follows from the standard result for energies expressed in
       * isochoric principal stretches (e.g. Holzapfel, "Nonlinear Solid Mechanics", Eq. 6.98):
       * \f[
       *   \tau_i = \sum_{p=1}^N \mu_p \bar\lambda_i^{\alpha_p} - \frac{1}{3} \sum_{j=1}^3
       *   \sum_{p=1}^N \mu_p \bar\lambda_j^{\alpha_p}\, ,
       * \f]
       * which is trace-free (deviatoric) by construction, consistent with the absence of a
       * volumetric term.
       *
       * @tparam T Scalar type, e.g. double, autodiff::dual.
       * @param principalStretches Principal stretches \f$\lambda_i\f$ of \f$\boldsymbol{F}\f$
       * (e.g. from DeformationMeasures::principalStretches()).
       * @param J Volume ratio \f$J = \det(\boldsymbol{F})\f$.
       * @param mu Array of \p nTerms Ogden moduli \f$\mu_p\f$.
       * @param alpha Array of \p nTerms Ogden exponents \f$\alpha_p\f$.
       * @param nTerms Number of Ogden terms \f$N\f$.
       * @return A pair containing the isochoric energy density and the principal Kirchhoff stress
       * \f$\tau_i\f$, diagonal in the same principal basis as \p principalStretches.
       */
      template < typename T >
      std::pair< T, Tensor3t< T > > OgdenPrincipalKirchhoffStress( const Tensor3t< T >& principalStretches,
                                                                   const T&             J,
                                                                   const double*        mu,
                                                                   const double*        alpha,
                                                                   int                  nTerms )
      {
        const T Jm13 = pow( J, -1. / 3. );

        Tensor3t< T > lambdaBar;
        for ( int i = 0; i < 3; ++i )
          lambdaBar( i ) = Jm13 * principalStretches( i );

        T             psi( 0. );
        Tensor3t< T > tauSum;
        for ( int i = 0; i < 3; ++i )
          tauSum( i ) = T( 0. );

        for ( int p = 0; p < nTerms; ++p ) {
          const double mu_p    = mu[p];
          const double alpha_p = alpha[p];
          T            sumPow( 0. );
          for ( int i = 0; i < 3; ++i ) {
            const T lambdaBarPow = pow( lambdaBar( i ), alpha_p );
            tauSum( i ) += mu_p * lambdaBarPow;
            sumPow += lambdaBarPow;
          }
          psi += ( mu_p / alpha_p ) * ( sumPow - 3. );
        }

        const T meanTau = ( tauSum( 0 ) + tauSum( 1 ) + tauSum( 2 ) ) / 3.;

        Tensor3t< T > tauPrincipal;
        for ( int i = 0; i < 3; ++i )
          tauPrincipal( i ) = tauSum( i ) - meanTau;

        return { psi, tauPrincipal };
      }

      /** @brief Classical isochoric Neo-Hookean hyperelastic potential.
       *
       * \f[
       *   \Psi_\mathrm{iso}(\boldsymbol C) = \frac{\mu}{2} \left( I_1 J^{-2/3} - 3 \right)
       * \f]
       * with \f$I_1=\operatorname{tr}\boldsymbol C\f$ and \f$J=\sqrt{\det\boldsymbol C}\f$. Unlike
       * Ogden's principal-stretch formulation, this potential is expressed directly in terms of
       * \f$\boldsymbol C\f$'s invariants (it is the special case of a quadratic, i.e. integer,
       * exponent), so no spectral decomposition of \f$\boldsymbol C\f$ is required or performed -
       * matching the efficiency of the compressible potentials in EnergyDensityFunctions::Compressible.
       *
       * @tparam T Scalar type, e.g. double, autodiff::dual.
       * @param C Right Cauchy-Green tensor.
       * @param mu Shear modulus.
       * @return Energy density.
       */
      template < typename T >
      T NeoHookePotential( const Tensor33t< T >& C, const double mu )
      {
        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );

        return mu / 2. * ( I1 * pow( J, -2. / 3. ) - 3. );
      }

      /** @brief Classical isochoric Mooney-Rivlin hyperelastic potential.
       *
       * \f[
       *   \Psi_\mathrm{iso}(\boldsymbol C) = C_1 \left( \bar I_1 - 3 \right) + C_2 \left( \bar I_2 -
       *   3 \right)
       * \f]
       * with \f$\bar I_1 = I_1 J^{-2/3}\f$, \f$\bar I_2 = I_2 J^{-4/3}\f$,
       * \f$I_1=\operatorname{tr}\boldsymbol C\f$, \f$I_2=\tfrac12\left(I_1^2-\operatorname{tr}
       * \boldsymbol C^2\right)\f$ and \f$J=\sqrt{\det\boldsymbol C}\f$. As for
       * NeoHookePotential(), no spectral decomposition of \f$\boldsymbol C\f$ is required.
       *
       * @tparam T Scalar type, e.g. double, autodiff::dual.
       * @param C Right Cauchy-Green tensor.
       * @param C1 First Mooney-Rivlin modulus.
       * @param C2 Second Mooney-Rivlin modulus.
       * @return Energy density.
       */
      template < typename T >
      T MooneyRivlinPotential( const Tensor33t< T >& C, const double C1, const double C2 )
      {
        const T J  = sqrt( determinant( C ) );
        const T I1 = trace( C );
        const T I2 = 0.5 * ( I1 * I1 - trace( C % C ) );

        return C1 * ( I1 * pow( J, -2. / 3. ) - 3. ) + C2 * ( I2 * pow( J, -4. / 3. ) - 3. );
      }

      namespace SecondOrderDerived {

        /** @brief Classical isochoric Neo-Hookean hyperelastic potential and its first and second
         * derivative w.r.t. \f$\boldsymbol C\f$.
         *
         * @details Mirrors EnergyDensityFunctions::Compressible::SecondOrderDerived::PenceGouPotentialB(),
         * dropping its volumetric (\f$K\f$-dependent) terms - the \f$J\f$/\f$I_1\f$ chain-rule
         * building blocks (\f$\partial J/\partial\boldsymbol C\f$,
         * \f$\partial^2 J/\partial\boldsymbol C\partial\boldsymbol C\f$,
         * \f$\partial I_1/\partial\boldsymbol C\f$) are identical.
         *
         * @tparam T Scalar type, e.g. double, float.
         * @param C Right Cauchy-Green tensor.
         * @param mu Shear modulus.
         * @return A tuple containing the energy density, its first and second derivative w.r.t. C.
         */
        template < typename T >
        std::tuple< T, Tensor33t< T >, Tensor3333t< T > > NeoHookePotential( const Tensor33t< T >& C, const double mu )
        {
          using namespace FastorIndices;

          const T J  = sqrt( determinant( C ) );
          const T I1 = trace( C );

          const T psi = mu / 2. * ( I1 * pow( J, -2. / 3. ) - 3. );

          const T dPsi_dJ  = -mu / 3. * I1 * pow( J, -5. / 3. );
          const T dPsi_dI1 = mu / 2. * pow( J, -2. / 3. );

          const Tensor33t< T > CInv   = inverse( C );
          const Tensor33t< T > dJ_dC  = 0.5 * J * transpose( CInv );
          const Tensor33t< T > dI1_dC = Spatial3D::I;

          const Tensor33t< T > dPsi_dC = dPsi_dJ * dJ_dC + dPsi_dI1 * dI1_dC;

          const T d2Psi_dJdJ  = 5. / 9. * mu * I1 * pow( J, -8. / 3. );
          const T d2Psi_dJdI1 = -mu / 3. * pow( J, -5. / 3. );

          const Tensor3333t< T > d2J_dCdC = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                            J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );

          const Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) + dPsi_dJ * d2J_dCdC +
                                              d2Psi_dJdI1 * ( einsum< IJ, KL >( dJ_dC, dI1_dC ) +
                                                              einsum< IJ, KL >( dI1_dC, dJ_dC ) );

          return { psi, dPsi_dC, d2Psi_dCdC };
        }

        /** @brief Classical isochoric Mooney-Rivlin hyperelastic potential and its first and second
         * derivative w.r.t. \f$\boldsymbol C\f$.
         *
         * @details Reuses the same \f$J\f$/\f$I_1\f$ building blocks as NeoHookePotential() (in
         * turn identical to those in
         * EnergyDensityFunctions::Compressible::SecondOrderDerived::PenceGouPotentialB()), adding
         * the second invariant \f$I_2=\tfrac12(I_1^2-\operatorname{tr}\boldsymbol C^2)\f$ with
         * \f$\partial I_2/\partial\boldsymbol C = I_1\boldsymbol I-\boldsymbol C\f$ and
         * \f$\partial^2 I_2/\partial\boldsymbol C\partial\boldsymbol C = \boldsymbol I\otimes
         * \boldsymbol I - \boldsymbol I^{\rm sym}\f$ (\f$\boldsymbol I^{\rm sym}\f$ the symmetric
         * fourth-order identity, i.e. \f$\partial\boldsymbol C/\partial\boldsymbol C\f$ for
         * symmetric \f$\boldsymbol C\f$).
         *
         * @tparam T Scalar type, e.g. double, float.
         * @param C Right Cauchy-Green tensor.
         * @param C1 First Mooney-Rivlin modulus.
         * @param C2 Second Mooney-Rivlin modulus.
         * @return A tuple containing the energy density, its first and second derivative w.r.t. C.
         */
        template < typename T >
        std::tuple< T, Tensor33t< T >, Tensor3333t< T > > MooneyRivlinPotential( const Tensor33t< T >& C,
                                                                                 const double          C1,
                                                                                 const double          C2 )
        {
          using namespace FastorIndices;

          const T J  = sqrt( determinant( C ) );
          const T I1 = trace( C );
          const T I2 = 0.5 * ( I1 * I1 - trace( C % C ) );

          const T psi = C1 * ( I1 * pow( J, -2. / 3. ) - 3. ) + C2 * ( I2 * pow( J, -4. / 3. ) - 3. );

          const T dPsi_dJ  = C1 * I1 * ( -2. / 3. ) * pow( J, -5. / 3. ) + C2 * I2 * ( -4. / 3. ) * pow( J, -7. / 3. );
          const T dPsi_dI1 = C1 * pow( J, -2. / 3. );
          const T dPsi_dI2 = C2 * pow( J, -4. / 3. );

          const Tensor33t< T > CInv   = inverse( C );
          const Tensor33t< T > dJ_dC  = 0.5 * J * transpose( CInv );
          const Tensor33t< T > dI1_dC = Spatial3D::I;
          const Tensor33t< T > dI2_dC = I1 * dI1_dC - C;

          const Tensor33t< T > dPsi_dC = dPsi_dJ * dJ_dC + dPsi_dI1 * dI1_dC + dPsi_dI2 * dI2_dC;

          const T d2Psi_dJdJ  = C1 * I1 * 10. / 9. * pow( J, -8. / 3. ) + C2 * I2 * 28. / 9. * pow( J, -10. / 3. );
          const T d2Psi_dJdI1 = C1 * ( -2. / 3. ) * pow( J, -5. / 3. );
          const T d2Psi_dJdI2 = C2 * ( -4. / 3. ) * pow( J, -7. / 3. );

          const Tensor3333t< T > d2J_dCdC = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                            J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );
          const Tensor3333t< T > d2I2_dCdC = einsum< IJ, KL >( dI1_dC, dI1_dC ) - Spatial3D::ISymm;

          const Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) + dPsi_dJ * d2J_dCdC +
                                              dPsi_dI2 * d2I2_dCdC +
                                              d2Psi_dJdI1 * ( einsum< IJ, KL >( dJ_dC, dI1_dC ) +
                                                              einsum< IJ, KL >( dI1_dC, dJ_dC ) ) +
                                              d2Psi_dJdI2 * ( einsum< IJ, KL >( dJ_dC, dI2_dC ) +
                                                              einsum< IJ, KL >( dI2_dC, dJ_dC ) );

          return { psi, dPsi_dC, d2Psi_dCdC };
        }

      } // namespace SecondOrderDerived

    }   // namespace Incompressible

  }     // namespace EnergyDensityFunctions

} // namespace Marmot::ContinuumMechanics
