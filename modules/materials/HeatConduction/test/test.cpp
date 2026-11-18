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
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTransportMaterialFactory.h"
#include <cmath>
#include <memory>
#include <stdexcept>

using namespace Marmot::Testing;
using Marmot::Materials::TransportMaterial;

namespace {

  std::unique_ptr< TransportMaterial > create( const double* p, int n )
  {
    return std::unique_ptr< TransportMaterial >(
      Marmot::Registry::TransportMaterialFactory::createMaterial( "HEATCONDUCTION", p, n, 1 ) );
  }

  void evaluate( const TransportMaterial&               material,
                 const TransportMaterial::FieldState&   fs,
                 TransportMaterial::Response&           r,
                 TransportMaterial::AlgorithmicTangent& t )
  {
    material.initializeResponse( r, t );
    r.stateVars = nullptr;
    material.computeTransport( r, t, fs, { 0.0, 1.0 } );
  }

  void testShortFormIsotropicConstant()
  {
    const double p[]      = { 2400., 900., 2.5 };
    auto         material = create( p, 3 );
    throwExceptionOnFailure( material->getNumberOfFields() == 1, "one field" );

    TransportMaterial::FieldState fs;
    fs.field    = TransportMaterial::FieldVector::Constant( 1, 30.0 );
    fs.fieldOld = TransportMaterial::FieldVector::Constant( 1, 20.0 );
    fs.gradient.setZero( 3 );
    fs.gradient << 1.0, -2.0, 3.0;
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    evaluate( *material, fs, r, t );

    throwExceptionOnFailure( checkIfEqual( r.contentIncrement( 0 ), 2400. * 900. * 10.0, 1e-8 ), "enthalpy increment" );
    throwExceptionOnFailure( checkIfEqual( r.flux( 0 ), -2.5, 1e-14 ) && checkIfEqual( r.flux( 1 ), 5.0, 1e-14 ) &&
                               checkIfEqual( r.flux( 2 ), -7.5, 1e-14 ),
                             "Fourier flux" );
  }

  // finite differences of all tangents of the full (anisotropic, temperature dependent) form
  void testTangentsByFiniteDifferences()
  {
    const double p[]      = { 2400., 900., 2.5, 1.5, 1.0, 0.3, -0.2, 0.1, 0.004, 0.002, 20.0 };
    auto         material = create( p, 11 );

    TransportMaterial::FieldState fs;
    fs.field    = TransportMaterial::FieldVector::Constant( 1, 65.0 );
    fs.fieldOld = TransportMaterial::FieldVector::Constant( 1, 40.0 );
    fs.gradient.setZero( 3 );
    fs.gradient << 120.0, -45.0, 80.0;

    TransportMaterial::Response           r, rp, rm;
    TransportMaterial::AlgorithmicTangent t, tp;
    evaluate( *material, fs, r, t );

    const double h = 1e-5;
    {
      auto fsp = fs, fsm = fs;
      fsp.field( 0 ) += h;
      fsm.field( 0 ) -= h;
      evaluate( *material, fsp, rp, tp );
      evaluate( *material, fsm, rm, tp );
      const double dH = ( rp.contentIncrement( 0 ) - rm.contentIncrement( 0 ) ) / ( 2 * h );
      throwExceptionOnFailure( checkIfEqual( dH, t.dContentIncrement_dField( 0, 0 ), 1e-6 * std::abs( dH ) ),
                               "dContent/dT" );
      for ( int i = 0; i < 3; i++ ) {
        const double dq = ( rp.flux( i ) - rm.flux( i ) ) / ( 2 * h );
        throwExceptionOnFailure( checkIfEqual( dq, t.dFlux_dField( i, 0 ), 1e-7 * ( 1 + std::abs( dq ) ) ),
                                 "dFlux/dT" );
      }
    }
    for ( int j = 0; j < 3; j++ ) {
      auto fsp = fs, fsm = fs;
      fsp.gradient( j ) += h;
      fsm.gradient( j ) -= h;
      evaluate( *material, fsp, rp, tp );
      evaluate( *material, fsm, rm, tp );
      for ( int i = 0; i < 3; i++ ) {
        const double dq = ( rp.flux( i ) - rm.flux( i ) ) / ( 2 * h );
        throwExceptionOnFailure( checkIfEqual( dq, t.dFlux_dGradient( i, j ), 1e-7 * ( 1 + std::abs( dq ) ) ),
                                 "dFlux/dGradT" );
      }
    }

    // the enthalpy increment is H(T) - H(T_old) with H = rho c0 (theta + beta_c theta^2 / 2)
    const auto enthalpy = []( double T ) { return 2400. * 900. * ( ( T - 20. ) + 0.001 * ( T - 20. ) * ( T - 20. ) ); };
    const double expected = enthalpy( 65. ) - enthalpy( 40. );
    throwExceptionOnFailure( checkIfEqual( r.contentIncrement( 0 ), expected, 1e-9 * expected ), "enthalpy form" );
  }

  void testInvalidInput()
  {
    const auto throws = []( const double* p, int n ) {
      try {
        create( p, n );
      }
      catch ( const std::invalid_argument& ) {
        return true;
      }
      return false;
    };
    const double fourProperties[] = { 1., 1., 1., 1. };
    throwExceptionOnFailure( throws( fourProperties, 4 ), "4 properties must be rejected" );
    const double negativeCapacity[] = { 1., -1., 1. };
    throwExceptionOnFailure( throws( negativeCapacity, 3 ), "negative capacity must be rejected" );
    const double indefinite[] = { 1., 1., 1., 1., 1., 2., 0., 0., 0., 0., 0. };
    throwExceptionOnFailure( throws( indefinite, 11 ), "indefinite conductivity must be rejected" );

    // a conductivity factor 1 + beta_k (T - T_ref) <= 0 is reported, not used
    const double                  p[]      = { 1., 1., 1., -0.01, 0., 0. };
    auto                          material = create( p, 6 );
    TransportMaterial::FieldState fs;
    fs.field    = TransportMaterial::FieldVector::Constant( 1, 150.0 );
    fs.fieldOld = fs.field;
    fs.gradient.setZero( 3 );
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    bool                                  threw = false;
    try {
      evaluate( *material, fs, r, t );
    }
    catch ( const std::runtime_error& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "non-positive conductivity must throw" );
  }

} // namespace

int main()
{
  testShortFormIsotropicConstant();
  testTangentsByFiniteDifferences();
  testInvalidInput();
  return 0;
}
