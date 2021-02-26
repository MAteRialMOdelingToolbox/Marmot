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
#include "Marmot/HaighWestergaard.h"
#include <utility>

namespace Marmot {
  namespace ContinuumMechanics::CommonConstitutiveModels {

    class MenetreyWillam {
    public:
      struct MenetreyWillamParameters {
        double Af, Bf, Cf, m, e;
      } param;

      enum class MenetreyWillamType { Mises, Rankine, DruckerPrager, MohrCoulomb };

      MenetreyWillam( const double              ft,
                      const MenetreyWillamType& type = MenetreyWillamType::Rankine,
                      const double              fc   = 0 );

      void setParameters( const double ft, const double fc, const MenetreyWillamType& type );

      double        polarRadius( const double& theta ) const;
      static double polarRadius( const double& theta, const double& e );

      std::pair< double, double >        dPolarRadius_dTheta( const double& theta ) const;
      static std::pair< double, double > dPolarRadius_dTheta( const double& theta, const double& e );

      double yieldFunction( const HaighWestergaard::HaighWestergaardCoordinates& hw, const double varEps = 0.0 ) const;

      std::tuple< double, double, double > dYieldFunction_dHaighWestergaard(
        const ContinuumMechanics::HaighWestergaard::HaighWestergaardCoordinates& hw,
        const double                                                             varEps = 0.0 ) const;

      /*
       * Inline functions
       */

      static inline double abaqusMohrCoulombPotentialVarEpsToMenetreyWillam( const double varEps, const double psi )
      {
        return varEps * 2 * std::sin( psi );
      }

      static inline double e( const double fc, const double ft ) { return ( fc + 2 * ft ) / ( 2 * fc + ft ); }

      static inline double c( const double fc, const double ft )
      {
        const double phi_ = phi( fc, ft );
        return ft * ( 1 + std::sin( phi_ ) ) / ( 2 + std::cos( phi_ ) );
      }

      static inline double phi( const double fc, const double ft ) { return std::asin( ( fc - ft ) / ( fc + ft ) ); }

      static inline double ft( const double c, const double phi )
      {
        return 2 * c * std::cos( phi ) / ( 1 + std::sin( phi ) );
      }

      static inline double fc( const double c, const double phi )
      {
        return 2 * c * std::cos( phi ) / ( 1 - std::sin( phi ) );
      }
    };

  } // namespace ContinuumMechanics::CommonConstitutiveModels
} // namespace Marmot
