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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"

namespace Marmot::ContinuumMechanics {

  namespace EnergyDensityFunctions {

    using namespace Fastor;
    using namespace FastorStandardTensors;

    /** @brief Hyperelastic Energy Density Function Wa acc. Pence & Gou (2015), Eq. (2.11)
     *
     *  The energy density function \f$W_a\f$ is given as
     *  \f[
     *    W_a = \frac{G}{2} (I_1 - 3) + \left(\frac{K}{2} - \frac{G}{3}\right) (J - 1)^2 - G \ln(J)
     *  \f]
     *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
     *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
     * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
     * modulus, respectively.
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
     * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
     * modulus, respectively.
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

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );

      T res = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

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
     * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and shear
     * modulus, respectively.
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
       * \boldsymbol{C}} + \frac{\partial W_b}{\partial I_1} \frac{\partial I_1}{\partial \boldsymbol{C}} \f] where \f[
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
       *  The energy density function \f$W_b\f$ is given as
       *  \f[
       *    W_b = \frac{K}{8} \left(J - \frac{1}{J}\right)^2 + \frac{G}{2} \left(I_1 J^{-\frac{2}{3}} - 3\right)
       *  \f]
       *  where \f$ I_1 = \text{tr}(\boldsymbol{C}) \f$ is the first invariant of the right Cauchy-Green tensor
       *  \f$ \boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} \f$, \f$ J = \sqrt{\det(\boldsymbol{C})} =
       * \det(\boldsymbol{F}) \f$ is the determinant of the deformation gradient, and \f$ K, G \f$ are the bulk and
       * shear modulus, respectively.
       *
       *  Additionally, the first derivative w.r.t. \f$\boldsymbol{C}\f$ is computed as
       *  \f[
       *    \frac{\partial W_b}{\partial \boldsymbol{C}} = \frac{\partial W_b}{\partial J} \frac{\partial J}{\partial
       * \boldsymbol{C}} + \frac{\partial W_b}{\partial I_1} \frac{\partial I_1}{\partial \boldsymbol{C}} \f] where \f[
       *    \frac{\partial J}{\partial \boldsymbol{C}} = \frac{1}{2} J \boldsymbol{C}^{-1}
       *  \f]
       *  and
       *  \f[
       *    \frac{\partial I_1}{\partial \boldsymbol{C}} = \boldsymbol{I}
       *  \f]
       * The second derivative w.r.t. \f$\boldsymbol{C}\f$ is computed as
       * \f[
       *   \frac{\partial^2 W_b}{\partial \boldsymbol{C} \partial \boldsymbol{C}} =
       *   \frac{\partial^2 W_b}{\partial J^2} \frac{\partial J}{\partial \boldsymbol{C}} \otimes \frac{\partial
       * J}{\partial \boldsymbol{C}} + \frac{\partial W_b}{\partial J} \frac{\partial^2 J}{\partial \boldsymbol{C}
       * \partial \boldsymbol{C}} + \frac{\partial^2 W_b}{\partial J \partial I_1} \left( \frac{\partial J}{\partial
       * \boldsymbol{C}} \otimes \frac{\partial I_1}{\partial \boldsymbol{C}} + \frac{\partial I_1}{\partial
       * \boldsymbol{C}} \otimes \frac{\partial J}{\partial \boldsymbol{C}} \right) \f] where \f[ \frac{\partial^2
       * J}{\partial \boldsymbol{C} \partial \boldsymbol{C}} = \frac{J}{4} \left( \boldsymbol{C}^{-1} \otimes
       * \boldsymbol{C}^{-1} - 2 \, \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1} \right) \f] and \f[ \boldsymbol{A}
       * \odot \boldsymbol{B} = A_{iK} B_{jK} \f] is the symmetric tensor product.
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
        const T          d2Psi_dJdJ  = K / 4. * ( 1. + 3. / ( J * J * J * J ) ) + 5. / 9. * G * I1 * pow( J, -8. / 3. );
        const T          d2Psi_dJdI1 = -G / 3. * pow( J, -5. / 3. );
        Tensor3333t< T > d2J_dCdC    = J / 4. * einsum< JI, LK, to_IJKL >( CInv, CInv ) -
                                    J / 2. * einsum< JK, LI, to_IJKL >( CInv, CInv );

        Tensor3333t< T > d2Psi_dCdC = d2Psi_dJdJ * einsum< IJ, KL >( dJ_dC, dJ_dC ) + dPsi_dJ * d2J_dCdC +
                                      d2Psi_dJdI1 *
                                        ( einsum< IJ, KL >( dJ_dC, dI1_dC ) + einsum< IJ, KL >( dI1_dC, dJ_dC ) );

        return { psi, dPsi_dC, d2Psi_dCdC };
      }

    } // namespace SecondOrderDerived

  }   // namespace EnergyDensityFunctions

} // namespace Marmot::ContinuumMechanics
