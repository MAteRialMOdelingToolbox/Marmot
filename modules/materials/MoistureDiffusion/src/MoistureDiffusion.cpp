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
#include "Marmot/MoistureDiffusion.h"
#include "Marmot/MarmotJournal.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  MoistureDiffusion::MoistureDiffusion( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : TransportMaterial( materialProperties, nMaterialProperties, materialNumber ),
      D1( 0.0 ),
      alpha0( 1.0 ),
      hc( 0.5 ),
      n( 1.0 ),
      isotherm( Linear ),
      K( 1.0 ),
      C( 0.0 ),
      k( 0.0 ),
      Vm( 0.0 ),
      dwInf( 0.0 ),
      tauS( 1.0 )
  {
    const double* p = materialProperties;

    const auto reject = [&]( const std::string& what ) {
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << ": " << what );
    };

    if ( nMaterialProperties != 4 && nMaterialProperties != 6 && nMaterialProperties != 8 && nMaterialProperties != 10 )
      reject( Marmot::Utilities::MakeString()
              << "4, 6, 8 or 10 material properties expected, but " << nMaterialProperties << " were given." );

    D1     = p[0];
    alpha0 = p[1];
    hc     = p[2];
    n      = p[3];

    if ( nMaterialProperties >= 6 ) {
      if ( p[4] == 0.0 )
        isotherm = Linear;
      else if ( p[4] == 1.0 )
        isotherm = BSB;
      else
        reject( "the sorption isotherm type (property 4) must be 0 (linear) or 1 (BSB)." );

      if ( isotherm == Linear )
        K = p[5];
      else {
        if ( nMaterialProperties < 8 )
          reject( "the BSB isotherm needs C, k and V_m (properties 5, 6, 7)." );
        C  = p[5];
        k  = p[6];
        Vm = p[7];
      }
    }
    if ( nMaterialProperties == 10 ) {
      dwInf = p[8];
      tauS  = p[9];
    }

    if ( !( D1 > 0.0 ) )
      reject( "D_1 must be positive." );
    if ( !( alpha0 > 0.0 && alpha0 <= 1.0 ) )
      reject( "alpha_0 must be in (0, 1]." );
    if ( !( hc > 0.0 && hc < 1.0 ) )
      reject( "h_c must be in (0, 1)." );
    if ( !( n >= 1.0 ) )
      reject( "n must be at least 1." );
    if ( isotherm == Linear && !( K >= 0.0 ) )
      reject( "the moisture capacity K must not be negative." );
    if ( isotherm == BSB && !( C > 0.0 && k > 0.0 && k < 1.0 && Vm > 0.0 ) )
      reject( "the BSB isotherm needs C > 0, 0 < k < 1 and V_m > 0." );
    if ( !( dwInf >= 0.0 ) )
      reject( "the self-desiccation water consumption must not be negative." );
    if ( dwInf > 0.0 && !( tauS > 0.0 ) )
      reject( "the self-desiccation time constant must be positive." );

    stateLayout.finalize();
  }

  double MoistureDiffusion::diffusivity( double h, double& dD_dh ) const
  {
    if ( h >= 1.0 ) {
      dD_dh = 0.0;
      return D1;
    }
    const double u   = ( 1.0 - h ) / ( 1.0 - hc );
    const double un1 = std::pow( u, n - 1.0 );
    const double g   = 1.0 + un1 * u;
    dD_dh            = D1 * ( 1.0 - alpha0 ) * n * un1 / ( ( 1.0 - hc ) * g * g );
    return D1 * ( alpha0 + ( 1.0 - alpha0 ) / g );
  }

  double MoistureDiffusion::waterContent( double h, double& dw_dh ) const
  {
    if ( isotherm == Linear ) {
      dw_dh = K;
      return K * h;
    }
    const double A = 1.0 - k * h;
    const double B = 1.0 + ( C - 1.0 ) * k * h;
    if ( !( A > 0.0 && B > 0.0 ) )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": the BSB isotherm is not defined at h = " << h << "." );
    const double c0    = C * k * Vm;
    const double f     = A * B;
    const double df_dh = -k * B + ( C - 1.0 ) * k * A;
    dw_dh              = c0 * ( f - h * df_dh ) / ( f * f );
    return c0 * h / f;
  }

  double MoistureDiffusion::selfDesiccationSource( double time, double dT ) const
  {
    if ( dwInf == 0.0 )
      return 0.0;
    if ( !( dT > 0.0 ) )
      return -dwInf / tauS * std::exp( -time / tauS );
    // exact average over [time - dT, time]: dwInf (exp(-(t - dT) / tau) - exp(-t / tau)) / dT
    return -dwInf * std::exp( -time / tauS ) * std::expm1( dT / tauS ) / dT;
  }

  void MoistureDiffusion::computeTransport( Response&            response,
                                            AlgorithmicTangent&  tangent,
                                            const FieldState&    fieldState,
                                            const TimeIncrement& timeIncrement ) const
  {
    const double h    = fieldState.field( 0 );
    const double hOld = fieldState.fieldOld( 0 );

    double       dw_dh = 0.0, dwOld_dh = 0.0, dD_dh = 0.0;
    const double w    = waterContent( h, dw_dh );
    const double wOld = waterContent( hOld, dwOld_dh );
    const double D    = diffusivity( h, dD_dh );

    const Eigen::Vector3d gradH = fieldState.gradient.head< 3 >();

    response.contentIncrement( 0 )                = w - wOld;
    tangent.dContentIncrement_dField( 0, 0 )      = dw_dh;
    response.flux.head< 3 >()                     = -D * gradH;
    tangent.dFlux_dField.block< 3, 1 >( 0, 0 )    = -dD_dh * gradH;
    tangent.dFlux_dGradient.block< 3, 3 >( 0, 0 ) = -D * Eigen::Matrix3d::Identity();
    response.source( 0 )                          = selfDesiccationSource( timeIncrement.time, timeIncrement.dT );
    tangent.dSource_dField( 0, 0 )                = 0.0;
  }

} // namespace Marmot::Materials
