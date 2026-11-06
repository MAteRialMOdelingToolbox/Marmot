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
#include "Marmot/ModelCode2010BondSlip.h"
#include "Marmot/MarmotJournal.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  namespace {
    double computeS0( const double* p, int nP )
    {
      if ( nP < 8 )
        throw std::invalid_argument( Marmot::Utilities::MakeString() << "ModelCode2010BondSlip: requires 8 properties: "
                                                                        "tauMax, s1, s2, s3, alpha, tauF, K0, Kn" );
      const double tauMax = p[0], s1 = p[1], alpha = p[4], K0 = p[6];
      if ( alpha >= 1.0 )
        return 0.0;
      return std::pow( tauMax / ( K0 * std::pow( s1, alpha ) ), 1. / ( 1. - alpha ) );
    }
  } // namespace

  ModelCode2010BondSlip::ModelCode2010BondSlip( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber )
    : Marmot::Materials::BondSlipLaw( materialProperties, nMaterialProperties, materialNumber ),
      tauMax( materialProperties[0] ),
      s1( materialProperties[1] ),
      s2( materialProperties[2] ),
      s3( materialProperties[3] ),
      alpha( materialProperties[4] ),
      tauF( materialProperties[5] ),
      K0( materialProperties[6] ),
      Kn( materialProperties[7] ),
      s0( computeS0( materialProperties, nMaterialProperties ) )
  {
    if ( !( s1 > 0 && s1 <= s2 && s2 < s3 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": requires 0 < s1 <= s2 < s3" );
    if ( !( alpha >= 0 && alpha <= 1 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": requires 0 <= alpha <= 1" );
    if ( !( tauF >= 0 && tauF <= tauMax ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": requires 0 <= tauF <= tauMax" );
    if ( !( K0 > tauMax / s1 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": requires K0 > tauMax / s1, the elastic branch must be steeper "
                                      "than the secant to the peak" );
    if ( !( Kn >= 0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": requires Kn >= 0" );

    stateLayout.add( "plastic slip", 1 );
    stateLayout.add( "max slip", 1 );
    stateLayout.finalize();
  }

  std::pair< double, double > ModelCode2010BondSlip::envelope( double s ) const
  {
    if ( s <= s0 )
      return { K0 * s, K0 };
    if ( s <= s1 ) {
      const double tau = tauMax * std::pow( s / s1, alpha );
      return { tau, alpha * tau / s };
    }
    if ( s <= s2 )
      return { tauMax, 0.0 };
    if ( s <= s3 ) {
      const double slope = -( tauMax - tauF ) / ( s3 - s2 );
      return { tauMax + slope * ( s - s2 ), slope };
    }
    return { tauF, 0.0 };
  }

  void ModelCode2010BondSlip::computeBondStress( State&                        state,
                                                 Marmot::Math::Matrix3d&       dBondStress_dSlip,
                                                 const Marmot::Math::Vector3d& slip,
                                                 const Marmot::Math::Vector3d&,
                                                 const TimeIncrement& ) const
  {
    double& sP    = stateLayout.getAs< double& >( state.stateVars, "plastic slip" );
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "max slip" );

    dBondStress_dSlip.setZero();

    // normal directions: linear elastic
    dBondStress_dSlip( 1, 1 ) = Kn;
    dBondStress_dSlip( 2, 2 ) = Kn;
    state.bondStress( 1 )     = Kn * slip( 1 );
    state.bondStress( 2 )     = Kn * slip( 2 );

    // tangential direction: elastic predictor ...
    const double s        = slip( 0 );
    const double tauOld   = state.bondStress( 0 );
    const double sPOld    = sP;
    const double kappaOld = kappa;

    kappa                          = std::max( kappaOld, std::abs( s ) );
    const auto [bound, boundSlope] = envelope( kappa );

    const double tauTrial = K0 * ( s - sPOld );

    if ( std::abs( tauTrial ) <= bound * ( 1 + 1e-12 ) ) {
      state.bondStress( 0 )     = tauTrial;
      dBondStress_dSlip( 0, 0 ) = K0;
    }
    else {
      // ... and frictional slip on the envelope bound
      const double direction = tauTrial > 0 ? 1.0 : -1.0;
      state.bondStress( 0 )  = direction * bound;
      sP                     = s - state.bondStress( 0 ) / K0;

      // the bound grows (or softens) with the slip only while the largest slip is exceeded
      const bool virginSlip     = std::abs( s ) > kappaOld;
      dBondStress_dSlip( 0, 0 ) = virginSlip ? direction * boundSlope * ( s > 0 ? 1.0 : -1.0 ) : 0.0;
    }

    state.dissipation += 0.5 * ( state.bondStress( 0 ) + tauOld ) * ( sP - sPOld );
    state.elasticEnergyDensity = 0.5 * state.bondStress( 0 ) * state.bondStress( 0 ) / K0 +
                                 0.5 * Kn * ( slip( 1 ) * slip( 1 ) + slip( 2 ) * slip( 2 ) );
  }
} // namespace Marmot::Materials
