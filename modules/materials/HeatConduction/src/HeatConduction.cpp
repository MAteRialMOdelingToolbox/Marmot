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
#include "Marmot/HeatConduction.h"
#include "Marmot/MarmotJournal.h"
#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  HeatConduction::HeatConduction( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : TransportMaterial( materialProperties, nMaterialProperties, materialNumber ),
      rhoC0( 0.0 ),
      k0( Marmot::Math::Matrix3d::Zero() ),
      betaK( 0.0 ),
      betaC( 0.0 ),
      TReference( 0.0 )
  {
    const double* p = materialProperties;

    switch ( nMaterialProperties ) {
    case 3:
    case 6: {
      k0 = p[2] * Marmot::Math::Matrix3d::Identity();
      if ( nMaterialProperties == 6 ) {
        betaK      = p[3];
        betaC      = p[4];
        TReference = p[5];
      }
      break;
    }
    case 11: {
      k0 << p[2], p[5], p[6], //
        p[5], p[3], p[7],     //
        p[6], p[7], p[4];
      betaK      = p[8];
      betaC      = p[9];
      TReference = p[10];
      break;
    }
    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": 3, 6 or 11 material properties expected, but "
                                   << nMaterialProperties << " were given." );
    }

    rhoC0 = p[0] * p[1];
    if ( !( rhoC0 >= 0.0 ) || !std::isfinite( rhoC0 ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__
                                                                   << ": the heat capacity rho * c must be finite "
                                                                      "and not negative." );

    if ( !std::isfinite( betaK ) || !std::isfinite( betaC ) || !std::isfinite( TReference ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": beta_k, beta_c and T_ref must be finite." );

    // the Cholesky decomposition does not report a NaN or infinite conductivity
    if ( !k0.allFinite() || Eigen::LLT< Marmot::Math::Matrix3d >( k0 ).info() != Eigen::Success )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": the conductivity must be finite and positive definite." );

    stateLayout.finalize();
  }

  double HeatConduction::enthalpy( double T ) const
  {
    const double theta = T - TReference;
    return rhoC0 * ( theta + 0.5 * betaC * theta * theta );
  }

  void HeatConduction::computeTransport( Response&                             response,
                                         AlgorithmicTangent&                   tangent,
                                         const FieldState&                     fieldState,
                                         [[maybe_unused]] const TimeIncrement& timeIncrement ) const
  {
    const double T     = fieldState.field( 0 );
    const double theta = T - TReference;

    const double conductivityFactor = 1.0 + betaK * theta;
    if ( !( conductivityFactor > 0.0 ) )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__ << ": the conductivity factor 1 + beta_k (T - T_ref) = "
                                << conductivityFactor << " is not positive at T = " << T << "." );

    const Eigen::Vector3d gradT = fieldState.gradient.head< 3 >();

    response.contentIncrement( 0 )                = enthalpy( T ) - enthalpy( fieldState.fieldOld( 0 ) );
    tangent.dContentIncrement_dField( 0, 0 )      = rhoC0 * ( 1.0 + betaC * theta );
    response.flux.head< 3 >()                     = -conductivityFactor * k0 * gradT;
    tangent.dFlux_dField.block< 3, 1 >( 0, 0 )    = -betaK * k0 * gradT;
    tangent.dFlux_dGradient.block< 3, 3 >( 0, 0 ) = -conductivityFactor * k0;
    response.source( 0 )                          = 0.0;
    tangent.dSource_dField( 0, 0 )                = 0.0;
  }

} // namespace Marmot::Materials
