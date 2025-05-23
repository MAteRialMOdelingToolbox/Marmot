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

    template < typename T >
    // \brief Hyperelastic Energy Density Function Wa acc. Pence & Gou (2015), Eq. (2.11)
    T PenceGouPotentialA( const Tensor33t< T >& C, const double K, const double G )
    {

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );

      T res = G / 2. * ( I1 - 3. ) + ( K / 2. - G / 3. ) * pow( J - 1, 2 ) - G * log( J );

      return res;
    }

    template < typename T >
    // \brief Hyperelastic Energy Density Function Wb acc. Pence & Gou (2015), Eq. (2.12)
    T PenceGouPotentialB( const Tensor33t< T >& C, const double K, const double G )
    {

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );

      T res = K / 8. * pow( J - 1. / J, 2. ) + G / 2. * ( I1 * pow( J, -2. / 3 ) - 3. );

      return res;
    }

    template < typename T >
    // \brief Hyperelastic Energy Density Function Wc acc. Pence & Gou (2015), Eq. (2.13)
    T PenceGouPotentialC( const Tensor33t< T >& C, const double K, const double G )
    {

      const T J  = sqrt( determinant( C ) );
      const T I1 = trace( C );

      T res = G / 2. * ( I1 - 3. ) + 3. * G * G / ( 3. * K - 2. * G ) * ( pow( J, 2. / 3 - K / G ) - 1 );

      return res;
    }

    namespace FirstOrderDerived {

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
