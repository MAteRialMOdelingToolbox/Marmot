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
#include "Marmot/ModelCode2010CoulombBondSlip.h"
#include "Marmot/MarmotJournal.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  namespace {
    double computeS0( const double* p, int nP )
    {
      if ( nP < 10 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "ModelCode2010CoulombBondSlip: requires 10 properties: "
                                        "tauMax, s1, s2, s3, alpha, tauF, K0, Kn, mu, pMax" );
      const double tauMax = p[0], s1 = p[1], alpha = p[4], K0 = p[6];
      if ( alpha >= 1.0 )
        return 0.0;
      return std::pow( tauMax / ( K0 * std::pow( s1, alpha ) ), 1. / ( 1. - alpha ) );
    }
  } // namespace

  ModelCode2010CoulombBondSlip::ModelCode2010CoulombBondSlip( const double* materialProperties,
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
      mu( materialProperties[8] ),
      pMax( materialProperties[9] ),
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
    if ( !( Kn > 0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": requires Kn > 0" );
    if ( !( mu >= 0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": requires mu >= 0" );

    stateLayout.add( "plastic slip", 1 );
    stateLayout.add( "max slip", 1 );
    stateLayout.add( "plastic normal slip", 2 );
    stateLayout.finalize();
  }

  std::pair< double, double > ModelCode2010CoulombBondSlip::envelope( double s ) const
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

  void ModelCode2010CoulombBondSlip::computeBondStress( State&                        state,
                                                        Marmot::Math::Matrix3d&       dBondStress_dSlip,
                                                        const Marmot::Math::Vector3d& slip,
                                                        const Marmot::Math::Vector3d&,
                                                        const TimeIncrement& ) const
  {
    double& sP    = stateLayout.getAs< double& >( state.stateVars, "plastic slip" );
    double& kappa = stateLayout.getAs< double& >( state.stateVars, "max slip" );
    double* sNP   = stateLayout.getStateView( state.stateVars, "plastic normal slip" ).stateLocation;

    dBondStress_dSlip.setZero();

    // normal directions: elastic-perfectly plastic, isotropic cap |t_n| <= pMax in the normal plane
    const Eigen::Vector2d sN( slip( 1 ), slip( 2 ) );
    const Eigen::Vector2d sNPOld( sNP[0], sNP[1] );
    const Eigen::Vector2d tNTrial = Kn * ( sN - sNPOld );
    const double          pTrial  = tNTrial.norm();
    Eigen::Vector2d       tN      = tNTrial;
    Eigen::Matrix2d       dtN_dsN = Kn * Eigen::Matrix2d::Identity();
    Eigen::Vector2d       dp_dsN  = Eigen::Vector2d::Zero();
    if ( pMax > 0 && pTrial > pMax ) {
      const Eigen::Vector2d e = tNTrial / pTrial;
      tN                      = pMax * e;
      dtN_dsN                 = Kn * pMax / pTrial * ( Eigen::Matrix2d::Identity() - e * e.transpose() );
      // dp/dsN = 0 on the cap
    }
    else if ( pTrial > 0 )
      dp_dsN = Kn * tNTrial / pTrial;
    const Eigen::Vector2d sNPNew = sN - tN / Kn;
    const double          p      = tN.norm();

    state.bondStress( 1 )                            = tN( 0 );
    state.bondStress( 2 )                            = tN( 1 );
    dBondStress_dSlip.template block< 2, 2 >( 1, 1 ) = dtN_dsN;

    // tangential direction: elastic predictor ...
    const double s        = slip( 0 );
    const double tauOld   = state.bondStress( 0 );
    const double sPOld    = sP;
    const double kappaOld = kappa;

    kappa                        = std::max( kappaOld, std::abs( s ) );
    const auto [bond, bondSlope] = envelope( kappa );
    const double bound           = bond + mu * p;
    const double tauTrial        = K0 * ( s - sPOld );

    if ( std::abs( tauTrial ) <= bound * ( 1 + 1e-12 ) ) {
      state.bondStress( 0 )     = tauTrial;
      dBondStress_dSlip( 0, 0 ) = K0;
    }
    else {
      // ... and frictional slip on the bound, which depends on the largest slip and on the normal pressure
      const double direction = tauTrial > 0 ? 1.0 : -1.0;
      state.bondStress( 0 )  = direction * bound;
      sP                     = s - state.bondStress( 0 ) / K0;

      const bool virginSlip     = std::abs( s ) > kappaOld;
      dBondStress_dSlip( 0, 0 ) = virginSlip ? direction * bondSlope * ( s > 0 ? 1.0 : -1.0 ) : 0.0;
      dBondStress_dSlip( 0, 1 ) = direction * mu * dp_dsN( 0 );
      dBondStress_dSlip( 0, 2 ) = direction * mu * dp_dsN( 1 );
    }

    state.dissipation += 0.5 * ( state.bondStress( 0 ) + tauOld ) * ( sP - sPOld ) + pMax * ( sNPNew - sNPOld ).norm();
    sNP[0]                     = sNPNew( 0 );
    sNP[1]                     = sNPNew( 1 );
    state.elasticEnergyDensity = 0.5 * state.bondStress( 0 ) * state.bondStress( 0 ) / K0 + 0.5 * tN.squaredNorm() / Kn;
  }
} // namespace Marmot::Materials
