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
#include "Marmot/ThermoHygro.h"
#include "Marmot/MarmotJournal.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  namespace {
    /// the number of properties of the uncoupled heat conduction law (full form)
    constexpr int nHeatProperties = 11;
    /// the number of properties of the uncoupled moisture diffusion law (full form)
    constexpr int nMoistureProperties = 10;

    int checkedCount( int nMaterialProperties )
    {
      if ( nMaterialProperties != 21 && nMaterialProperties != 32 )
        throw std::invalid_argument( Marmot::Utilities::MakeString()
                                     << __PRETTY_FUNCTION__ << ": 21 or 32 material properties expected, but "
                                     << nMaterialProperties << " were given." );
      return nMaterialProperties;
    }
  } // namespace

  ThermoHygro::ThermoHygro( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : TransportMaterial( materialProperties, checkedCount( nMaterialProperties ), materialNumber ),
      heat( materialProperties, nHeatProperties, materialNumber ),
      moisture( materialProperties + nHeatProperties, nMoistureProperties, materialNumber ),
      TReference( materialProperties[10] ),
      dwInf( materialProperties[19] ),
      tauS( materialProperties[20] ),
      TZero( -273.15 ),
      QD( 0.0 ),
      mD( 0.0 ),
      kappa( 0.0 ),
      chi( 0.0 ),
      cw( 0.0 ),
      betaKh( 0.0 ),
      ell( 0.0 ),
      Eh( 0.0 ),
      QInf( 0.0 ),
      tauQ( 1.0 ),
      withMaturity( false )
  {
    const double* p = materialProperties;

    const auto reject = [&]( const std::string& what ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": " << what );
    };

    if ( nMaterialProperties == 32 ) {
      TZero  = p[21];
      QD     = p[22];
      mD     = p[23];
      kappa  = p[24];
      chi    = p[25];
      cw     = p[26];
      betaKh = p[27];
      ell    = p[28];
      Eh     = p[29];
      QInf   = p[30];
      tauQ   = p[31];
    }
    withMaturity = Eh != 0.0;

    if ( !( QD >= 0.0 ) )
      reject( "the activation temperature Q_D / R of the diffusivity must not be negative." );
    if ( !( chi >= 0.0 ) )
      reject( "the vapour pressure factor chi must not be negative." );
    if ( !( cw >= 0.0 ) )
      reject( "the specific heat of the water c_w must not be negative." );
    if ( !( betaKh > -1.0 ) )
      reject( "the moisture coefficient of the conductivity must be larger than -1." );
    if ( !( Eh >= 0.0 ) )
      reject( "the activation temperature E_h / R of the hydration must not be negative." );
    if ( !( QInf >= 0.0 ) )
      reject( "the heat of hydration must not be negative." );
    if ( QInf > 0.0 && !( tauQ > 0.0 ) )
      reject( "the time constant of the heat of hydration must be positive." );
    if ( !std::isfinite( kappa ) || !std::isfinite( mD ) || !std::isfinite( ell ) || !std::isfinite( TZero ) )
      reject( "kappa, m, ell and T_0 must be finite." );
    if ( ( QD != 0.0 || mD != 0.0 || chi != 0.0 || Eh != 0.0 ) && !( TReference - TZero > 0.0 ) )
      reject( "the absolute reference temperature T_ref - T_0 must be positive." );

    stateLayout.add( "equivalent age", 1 );
    stateLayout.finalize();
  }

  bool ThermoHygro::hasStorage( int field ) const
  {
    if ( field == 0 )
      return heat.hasStorage( 0 ) || ( cw > 0.0 && moisture.hasStorage( 0 ) );
    return moisture.hasStorage( 0 );
  }

  double ThermoHygro::absoluteTemperature( double T ) const
  {
    const double Ta = T - TZero;
    if ( !( Ta > 0.0 ) )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": non-positive absolute temperature at T = " << T << "." );
    return Ta;
  }

  double ThermoHygro::diffusivityTemperatureFactor( double T, double& df_dT ) const
  {
    if ( QD == 0.0 && mD == 0.0 ) {
      df_dT = 0.0;
      return 1.0;
    }
    const double Ta = absoluteTemperature( T ), TaRef = TReference - TZero;
    const double f = std::pow( Ta / TaRef, mD ) * std::exp( QD * ( 1.0 / TaRef - 1.0 / Ta ) );
    df_dT          = f * ( mD / Ta + QD / ( Ta * Ta ) );
    return f;
  }

  double ThermoHygro::maturityRate( double T, double& da_dT ) const
  {
    if ( !withMaturity ) {
      da_dT = 0.0;
      return 1.0;
    }
    const double Ta = absoluteTemperature( T ), TaRef = TReference - TZero;
    const double a = std::exp( Eh * ( 1.0 / TaRef - 1.0 / Ta ) );
    da_dT          = a * Eh / ( Ta * Ta );
    return a;
  }

  void ThermoHygro::computeTransport( Response&            response,
                                      AlgorithmicTangent&  tangent,
                                      const FieldState&    fieldState,
                                      const TimeIncrement& timeIncrement ) const
  {
    const double T = fieldState.field( 0 ), TOld = fieldState.fieldOld( 0 );
    const double h = fieldState.field( 1 ), hOld = fieldState.fieldOld( 1 );
    const double dT = timeIncrement.dT;

    const Eigen::Vector3d gradT = fieldState.gradient.segment< 3 >( 0 );
    const Eigen::Vector3d gradH = fieldState.gradient.segment< 3 >( 3 );
    const Eigen::Matrix3d I     = Eigen::Matrix3d::Identity();

    // ---- the uncoupled heat conduction law: dry enthalpy increment, Fourier flux
    Response           heatResponse;
    AlgorithmicTangent heatTangent;
    FieldState         heatState;
    heat.initializeResponse( heatResponse, heatTangent );
    heatState.field        = fieldState.field.head( 1 );
    heatState.fieldOld     = fieldState.fieldOld.head( 1 );
    heatState.gradient     = fieldState.gradient.head( 3 );
    heatResponse.stateVars = nullptr;
    heat.computeTransport( heatResponse, heatTangent, heatState, timeIncrement );

    // ---- water content: temperature dependent isotherm w(h - kappa (T - T_ref))
    double       dwIso = 0.0, dwIsoOld = 0.0;
    const double w     = moisture.waterContent( h - kappa * ( T - TReference ), dwIso );
    const double wOld  = moisture.waterContent( hOld - kappa * ( TOld - TReference ), dwIsoOld );
    const double dw_dh = dwIso, dw_dT = -kappa * dwIso;

    // ---- moisture flux q_w = -D_BN(h) f(T) (grad h + chi h g(T) grad T)
    double       dDBN_dh = 0.0, df_dT = 0.0;
    const double DBN = moisture.diffusivity( h, dDBN_dh );
    const double f   = diffusivityTemperatureFactor( T, df_dT );
    const double D   = DBN * f;

    double g = 0.0, dg_dT = 0.0;
    if ( chi != 0.0 ) {
      const double Ta = absoluteTemperature( T );
      g               = vapourTemperature / ( Ta * Ta );
      dg_dT           = -2.0 * g / Ta;
    }
    const double          chiG    = chi * h * g; // factor of grad T in the driving gradient
    const Eigen::Vector3d driving = gradH + chiG * gradT;

    const Eigen::Vector3d qW         = -D * driving;
    const Eigen::Vector3d dqW_dh     = -( dDBN_dh * f ) * driving - ( D * chi * g ) * gradT;
    const Eigen::Vector3d dqW_dT     = -( DBN * df_dT ) * driving - ( D * chi * h * dg_dT ) * gradT;
    const Eigen::Matrix3d dqW_dGradH = -D * I;
    const Eigen::Matrix3d dqW_dGradT = -( D * chiG ) * I;

    // ---- heat flux q_T = (1 + beta_kh h) q_Fourier(T) + ell q_w
    const double psi = 1.0 + betaKh * h;
    if ( !( psi > 0.0 ) )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": the moisture factor of the conductivity 1 + beta_kh h = "
                                << psi << " is not positive at h = " << h << "." );
    const Eigen::Vector3d qFourier = heatResponse.flux.head< 3 >();

    // ---- sources: heat of hydration and self-desiccation over the equivalent age
    double&      te    = stateLayout.getAs< double& >( response.stateVars, "equivalent age" );
    const double teN   = te;
    double       da_dT = 0.0;
    const double a     = maturityRate( T, da_dT );
    double       rT = 0.0, drT_dT = 0.0, rH = 0.0, drH_dT = 0.0;
    if ( !withMaturity ) {
      // the equivalent age is the analysis time: the same expressions as MoistureDiffusion
      te = timeIncrement.time;
      rH = moisture.selfDesiccationSource( timeIncrement.time, dT );
      if ( QInf > 0.0 )
        rT = dT > 0.0 ? QInf * std::exp( -timeIncrement.time / tauQ ) * std::expm1( dT / tauQ ) / dT
                      : QInf / tauQ * std::exp( -timeIncrement.time / tauQ );
    }
    else {
      const double dTe = a * dT;
      te               = teN + dTe;
      // exact average of A (1 - exp(-t_e / tau)) over the increment, and its derivative with respect to T
      const auto averageRate = [&]( double A, double tau, double& dRate_dT ) {
        if ( A == 0.0 ) {
          dRate_dT = 0.0;
          return 0.0;
        }
        const double e = std::exp( -te / tau );
        dRate_dT       = A * e / tau * da_dT;
        return dT > 0.0 ? A * e * std::expm1( dTe / tau ) / dT : A * e / tau * a;
      };
      rT              = averageRate( QInf, tauQ, drT_dT );
      double dSink_dT = 0.0;
      rH              = -averageRate( dwInf, tauS, dSink_dT );
      drH_dT          = -dSink_dT;
    }

    // ---- assemble: field 0 = T, field 1 = h
    const double dTIncrement = T - TOld;

    response.contentIncrement( 0 )           = heatResponse.contentIncrement( 0 ) + cw * w * dTIncrement;
    response.contentIncrement( 1 )           = w - wOld;
    tangent.dContentIncrement_dField( 0, 0 ) = heatTangent.dContentIncrement_dField( 0, 0 ) +
                                               cw * ( w + dw_dT * dTIncrement );
    tangent.dContentIncrement_dField( 0, 1 ) = cw * dw_dh * dTIncrement;
    tangent.dContentIncrement_dField( 1, 0 ) = dw_dT;
    tangent.dContentIncrement_dField( 1, 1 ) = dw_dh;

    response.flux.segment< 3 >( 0 ) = psi * qFourier + ell * qW;
    response.flux.segment< 3 >( 3 ) = qW;

    tangent.dFlux_dField.block< 3, 1 >( 0, 0 ) = psi * heatTangent.dFlux_dField.block< 3, 1 >( 0, 0 ) + ell * dqW_dT;
    tangent.dFlux_dField.block< 3, 1 >( 0, 1 ) = betaKh * qFourier + ell * dqW_dh;
    tangent.dFlux_dField.block< 3, 1 >( 3, 0 ) = dqW_dT;
    tangent.dFlux_dField.block< 3, 1 >( 3, 1 ) = dqW_dh;

    tangent.dFlux_dGradient.block< 3, 3 >( 0, 0 ) = psi * heatTangent.dFlux_dGradient.block< 3, 3 >( 0, 0 ) +
                                                    ell * dqW_dGradT;
    tangent.dFlux_dGradient.block< 3, 3 >( 0, 3 ) = ell * dqW_dGradH;
    tangent.dFlux_dGradient.block< 3, 3 >( 3, 0 ) = dqW_dGradT;
    tangent.dFlux_dGradient.block< 3, 3 >( 3, 3 ) = dqW_dGradH;

    response.source( 0 )           = rT;
    response.source( 1 )           = rH;
    tangent.dSource_dField( 0, 0 ) = drT_dT;
    tangent.dSource_dField( 0, 1 ) = 0.0;
    tangent.dSource_dField( 1, 0 ) = drH_dT;
    tangent.dSource_dField( 1, 1 ) = 0.0;
  }

} // namespace Marmot::Materials
