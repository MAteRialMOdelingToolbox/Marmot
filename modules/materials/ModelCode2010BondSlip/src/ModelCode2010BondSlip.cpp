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
    /// the property array, after checking the number of properties (before any property is read)
    const double* checkedProperties( const double* materialProperties, int nMaterialProperties )
    {
      if ( nMaterialProperties < 8 || nMaterialProperties > 10 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << "ModelCode2010BondSlip: requires 8 to 10 properties: "
                                        "tauMax, s1, s2, s3, alpha, tauF, K0, Kn, (mu, (pMax)), got "
                                     << nMaterialProperties );
      return materialProperties;
    }

    double optionalProperty( const double* materialProperties, int nMaterialProperties, int index )
    {
      return index < nMaterialProperties ? materialProperties[index] : 0.0;
    }

    /// where the linear regularization K0 s meets the power law tauMax (s / s1)^alpha
    double regularizationSlip( double tauMax, double s1, double alpha, double K0 )
    {
      if ( alpha >= 1.0 )
        return 0.0;
      return std::pow( tauMax / ( K0 * std::pow( s1, alpha ) ), 1. / ( 1. - alpha ) );
    }
  } // namespace

  ModelCode2010BondSlip::ModelCode2010BondSlip( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber )
    : Marmot::Materials::BondSlipLaw( materialProperties, nMaterialProperties, materialNumber ),
      tauMax( checkedProperties( materialProperties, nMaterialProperties )[0] ),
      s1( materialProperties[1] ),
      s2( materialProperties[2] ),
      s3( materialProperties[3] ),
      alpha( materialProperties[4] ),
      tauF( materialProperties[5] ),
      K0( materialProperties[6] ),
      Kn( materialProperties[7] ),
      mu( optionalProperty( materialProperties, nMaterialProperties, 8 ) ),
      pMax( optionalProperty( materialProperties, nMaterialProperties, 9 ) ),
      s0( regularizationSlip( tauMax, s1, alpha, K0 ) )
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
    if ( !( mu >= 0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": requires mu >= 0" );
    if ( pMax > 0 && !( Kn > 0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": a normal stress cap pMax > 0 requires Kn > 0" );

    stateLayout.add( "plastic slip", 1 );
    stateLayout.add( "max slip", 1 );
    stateLayout.add( "plastic normal slip", 2 );
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
                                                 const Marmot::Math::Vector3d& /*dSlip: rate independent*/,
                                                 const TimeIncrement& /*timeIncrement: rate independent*/ ) const
  {
    double&                       plasticSlip = stateLayout.getAs< double& >( state.stateVars, "plastic slip" );
    double&                       kappa       = stateLayout.getAs< double& >( state.stateVars, "max slip" );
    Eigen::Map< Eigen::Vector2d > plasticNormalSlip(
      stateLayout.getStateView( state.stateVars, "plastic normal slip" ).stateLocation );

    dBondStress_dSlip.setZero();

    // normal directions: elastic, optionally perfectly plastic with the isotropic cap |t_n| <= pMax
    const Eigen::Vector2d normalSlip( slip( 1 ), slip( 2 ) );
    const Eigen::Vector2d plasticNormalSlipOld = plasticNormalSlip;
    const Eigen::Vector2d normalStressTrial    = Kn * ( normalSlip - plasticNormalSlipOld );
    const double          pressureTrial        = normalStressTrial.norm();

    Eigen::Vector2d normalStress              = normalStressTrial;
    Eigen::Matrix2d dNormalStress_dNormalSlip = Kn * Eigen::Matrix2d::Identity();
    Eigen::Vector2d dPressure_dNormalSlip     = Eigen::Vector2d::Zero();
    double          crushingDissipation       = 0.0;
    const bool      normalCapActive           = pMax > 0 && pressureTrial > pMax;
    if ( normalCapActive ) {
      const Eigen::Vector2d normalDirection = normalStressTrial / pressureTrial;
      normalStress                          = pMax * normalDirection;
      dNormalStress_dNormalSlip             = Kn * pMax / pressureTrial *
                                  ( Eigen::Matrix2d::Identity() - normalDirection * normalDirection.transpose() );
      // dPressure_dNormalSlip = 0 on the cap
      plasticNormalSlip   = normalSlip - normalStress / Kn;
      crushingDissipation = pMax * ( plasticNormalSlip - plasticNormalSlipOld ).norm();
    }
    else if ( pressureTrial > 0 )
      dPressure_dNormalSlip = Kn * normalStressTrial / pressureTrial;
    const double pressure = normalStress.norm();

    state.bondStress( 1 )                            = normalStress( 0 );
    state.bondStress( 2 )                            = normalStress( 1 );
    dBondStress_dSlip.template block< 2, 2 >( 1, 1 ) = dNormalStress_dNormalSlip;

    // tangential direction: 1D elastoplastic return mapping
    //   yield function  f = |tau_trial| - tau_y(kappa_{n+1}, p),  tau_y = tau_env(kappa_{n+1}) + mu p
    //   hardening       kappa_{n+1} = max(kappa_n, |s_{n+1}|)
    // kappa is the largest *total* slip (not an accumulated plastic slip), so kappa_{n+1}, and with it the yield
    // stress, is known before the check: the return is closed form, no local Newton iteration is needed.
    //   f <= 0: tau = tau_trial = K0 (s - s_p),
    //   f > 0:  tau = sign(tau_trial) tau_y,  s_p = s - tau / K0.
    const double s              = slip( 0 );
    const double tauOld         = state.bondStress( 0 );
    const double plasticSlipOld = plasticSlip;
    const double kappaOld       = kappa;

    kappa                                         = std::max( kappaOld, std::abs( s ) );
    const auto [envelopeStress, dEnvelope_dKappa] = envelope( kappa );
    const double yieldBondStress                  = envelopeStress + mu * pressure;
    const double tauTrial                         = K0 * ( s - plasticSlipOld );
    const double yieldFunction                    = std::abs( tauTrial ) - yieldBondStress;

    if ( yieldFunction <= 1e-12 * yieldBondStress ) {
      state.bondStress( 0 )     = tauTrial;
      dBondStress_dSlip( 0, 0 ) = K0;
    }
    else {
      const double signTauTrial = tauTrial > 0 ? 1.0 : -1.0;
      state.bondStress( 0 )     = signTauTrial * yieldBondStress;
      plasticSlip               = s - state.bondStress( 0 ) / K0;

      // the envelope moves with the slip only while the largest slip is exceeded
      const double dKappa_dSlip = std::abs( s ) > kappaOld ? ( s > 0 ? 1.0 : -1.0 ) : 0.0;
      dBondStress_dSlip( 0, 0 ) = signTauTrial * dEnvelope_dKappa * dKappa_dSlip;
      // Coulomb friction: d(mu p)/ds_n
      dBondStress_dSlip( 0, 1 ) = signTauTrial * mu * dPressure_dNormalSlip( 0 );
      dBondStress_dSlip( 0, 2 ) = signTauTrial * mu * dPressure_dNormalSlip( 1 );
    }

    state.dissipation += 0.5 * ( state.bondStress( 0 ) + tauOld ) * ( plasticSlip - plasticSlipOld ) +
                         crushingDissipation;
    state.elasticEnergyDensity = 0.5 * state.bondStress( 0 ) * state.bondStress( 0 ) / K0 +
                                 0.5 * Kn * ( normalSlip - plasticNormalSlip ).squaredNorm();
  }
} // namespace Marmot::Materials
