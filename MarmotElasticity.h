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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Thomas Mader thomas.mader@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotVoigt.h"

namespace Marmot {

  namespace ContinuumMechanics {
    namespace Elasticity::Isotropic {

      double constexpr E( const double K, const double G ) { return 9. * K * G / ( 3. * K + G ); }

      double constexpr nu( const double K, const double G ) { return ( 3 * K - 2 * G ) / ( 6 * K + 2 * G ); }

      double constexpr shearModulus( const double E, const double nu ) { return E / ( 2 * ( 1 + nu ) ); }

      double constexpr lameParameter( const double E, const double nu )
      {
        return E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
      }

      Matrix6d complianceTensor( const double E, const double nu );
      Matrix6d stiffnessTensor( const double E, const double nu );

    } // namespace Elasticity::Isotropic

    namespace Elasticity::TransverseIsotropic {
      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double nu12,
                                 const double nu23,
                                 const double G12 );
      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double nu12,
                                const double nu23,
                                const double G12 );
    } // namespace Elasticity::TransverseIsotropic

    namespace Elasticity::Orthotropic {
      Matrix6d complianceTensor( const double E1,
                                 const double E2,
                                 const double E3,
                                 const double nu12,
                                 const double nu23,
                                 const double nu13,
                                 const double G12,
                                 const double G23,
                                 const double G31 );
      Matrix6d stiffnessTensor( const double E1,
                                const double E2,
                                const double E3,
                                const double nu12,
                                const double nu23,
                                const double nu13,
                                const double G12,
                                const double G23,
                                const double G31 );

    } // namespace Elasticity::Orthotropic
  }   // namespace ContinuumMechanics
} // namespace Marmot
