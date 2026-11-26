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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTransportMaterialFactory.h"
#include "Marmot/MoistureDiffusion.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace Marmot::Testing;
using Marmot::Materials::TransportMaterial;

namespace {

  std::unique_ptr< TransportMaterial > create( const double* p, int n )
  {
    return std::unique_ptr< TransportMaterial >(
      Marmot::Registry::TransportMaterialFactory::createMaterial( "MOISTUREDIFFUSION", p, n, 1 ) );
  }

  TransportMaterial::FieldState state( double h, double hOld, double gx, double gy, double gz )
  {
    TransportMaterial::FieldState fs;
    fs.field    = TransportMaterial::FieldVector::Constant( 1, h );
    fs.fieldOld = TransportMaterial::FieldVector::Constant( 1, hOld );
    fs.gradient.setZero( 3 );
    fs.gradient << gx, gy, gz;
    return fs;
  }

  void evaluate( const TransportMaterial&               material,
                 const TransportMaterial::FieldState&   fs,
                 TransportMaterial::Response&           r,
                 TransportMaterial::AlgorithmicTangent& t,
                 double                                 time = 10.0,
                 double                                 dT   = 2.0 )
  {
    material.initializeResponse( r, t );
    r.stateVars = nullptr;
    material.computeTransport( r, t, fs, { time, dT } );
  }

  // D1, alpha0, hc, n, BSB: C, k, Vm, self-desiccation dwInf, tauS
  const double fullBSB[] = { 2.0e-3, 0.05, 0.8, 15.0, 1.0, 18.0, 0.73, 20.0, 5.0, 30.0 };

  /// Bazant-Najjar diffusivity and linear isotherm in closed form
  void testClosedForms()
  {
    const double p[]      = { 3.0e-9, 0.05, 0.75, 6.0 };
    auto         material = create( p, 4 );
    throwExceptionOnFailure( material->getNumberOfFields() == 1 && material->hasStorage( 0 ), "one field, storage" );

    const double                          h = 0.6, hOld = 0.9;
    auto                                  fs = state( h, hOld, 1.0, -2.0, 0.5 );
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    evaluate( *material, fs, r, t );

    const double D = 3.0e-9 * ( 0.05 + 0.95 / ( 1.0 + std::pow( 0.4 / 0.25, 6.0 ) ) );
    throwExceptionOnFailure( checkIfEqual( r.contentIncrement( 0 ), h - hOld, 1e-15 ), "linear isotherm, K = 1" );
    throwExceptionOnFailure( checkIfEqual( r.flux( 0 ), -D, 1e-12 * D ) &&
                               checkIfEqual( r.flux( 1 ), 2.0 * D, 1e-12 * D ),
                             "Bazant-Najjar flux" );
    throwExceptionOnFailure( r.source( 0 ) == 0.0, "no sink" );

    // D(h_c) = D1 (1 + alpha0) / 2, D(1) = D1, D(h >= 1) = D1
    auto Dof = [&]( double hh ) {
      auto fs_ = state( hh, hh, 1.0, 0.0, 0.0 );
      evaluate( *material, fs_, r, t );
      return -r.flux( 0 );
    };
    throwExceptionOnFailure( checkIfEqual( Dof( 0.75 ), 3.0e-9 * 0.525, 1e-21 ), "D(h_c)" );
    throwExceptionOnFailure( checkIfEqual( Dof( 1.0 ), 3.0e-9, 1e-21 ), "D(1)" );
    throwExceptionOnFailure( checkIfEqual( Dof( 1.05 ), 3.0e-9, 1e-21 ), "D(h > 1)" );

    // linear isotherm without capacity: steady state
    const double steady[] = { 3.0e-9, 0.05, 0.75, 6.0, 0.0, 0.0 };
    throwExceptionOnFailure( !create( steady, 6 )->hasStorage( 0 ), "K = 0: no storage" );
  }

  /// BSB isotherm: content against the formula and the capacity by finite differences
  void testBSB()
  {
    auto       material = create( fullBSB, 10 );
    const auto w        = []( double h ) {
      return 18.0 * 0.73 * 20.0 * h / ( ( 1.0 - 0.73 * h ) * ( 1.0 + 17.0 * 0.73 * h ) );
    };
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    for ( double h : { 0.05, 0.3, 0.6, 0.85, 0.99, 1.0, 1.02 } ) {
      auto fs = state( h, 0.5, 0.0, 0.0, 0.0 );
      evaluate( *material, fs, r, t );
      throwExceptionOnFailure( checkIfEqual( r.contentIncrement( 0 ), w( h ) - w( 0.5 ), 1e-12 * w( 1.0 ) ),
                               "BSB content" );
      const double e    = 1e-6;
      const double dwFD = ( w( h + e ) - w( h - e ) ) / ( 2 * e );
      throwExceptionOnFailure( checkIfEqual( t.dContentIncrement_dField( 0, 0 ), dwFD, 1e-7 * dwFD ), "BSB capacity" );
    }
    // outside the range of the isotherm (k h >= 1): reported, not used
    bool threw = false;
    try {
      auto fs = state( 1.0 / 0.73 + 0.01, 0.5, 0.0, 0.0, 0.0 );
      evaluate( *material, fs, r, t );
    }
    catch ( const std::runtime_error& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "BSB beyond k h = 1 must throw" );
  }

  /// all tangents by central differences, over the range of h (around h_c, near 1, above 1), linear and BSB
  void testTangentsByFiniteDifferences()
  {
    const double linear[] = { 2.0e-3, 0.05, 0.8, 15.0, 0.0, 120.0 };
    const double nOne[]   = { 2.0e-3, 0.1, 0.7, 1.0 };
    double       maxErr   = 0.0;
    for ( auto [p, np] : { std::pair{ fullBSB, 10 }, std::pair{ linear, 6 }, std::pair{ nOne, 4 } } ) {
      auto material = create( p, np );
      for ( double h : { 0.1, 0.5, 0.75, 0.8, 0.85, 0.95, 0.999, 1.03 } ) {
        auto                                  fs = state( h, 0.62, 3.0, -1.5, 0.7 );
        TransportMaterial::Response           r, rp, rm;
        TransportMaterial::AlgorithmicTangent t, td;
        evaluate( *material, fs, r, t );

        const double e   = 1e-7;
        auto         fsp = fs, fsm = fs;
        fsp.field( 0 ) += e;
        fsm.field( 0 ) -= e;
        evaluate( *material, fsp, rp, td );
        evaluate( *material, fsm, rm, td );
        const double dw = ( rp.contentIncrement( 0 ) - rm.contentIncrement( 0 ) ) / ( 2 * e );
        const double ew = std::abs( dw - t.dContentIncrement_dField( 0, 0 ) ) / ( 1e-12 + std::abs( dw ) );
        double       eq = 0.0;
        for ( int i = 0; i < 3; i++ ) {
          const double dq = ( rp.flux( i ) - rm.flux( i ) ) / ( 2 * e );
          eq = std::max( eq, std::abs( dq - t.dFlux_dField( i, 0 ) ) / ( 1e-3 * 3.0 * 2e-3 + std::abs( dq ) ) );
        }
        const double es = std::abs( ( rp.source( 0 ) - rm.source( 0 ) ) / ( 2 * e ) - t.dSource_dField( 0, 0 ) );
        double       eg = 0.0;
        for ( int j = 0; j < 3; j++ ) {
          auto fgp = fs, fgm = fs;
          fgp.gradient( j ) += e;
          fgm.gradient( j ) -= e;
          evaluate( *material, fgp, rp, td );
          evaluate( *material, fgm, rm, td );
          for ( int i = 0; i < 3; i++ ) {
            const double dq = ( rp.flux( i ) - rm.flux( i ) ) / ( 2 * e );
            eg              = std::max( eg, std::abs( dq - t.dFlux_dGradient( i, j ) ) / 2e-3 );
          }
        }
        maxErr = std::max( { maxErr, ew, eq, eg, es } );
        throwExceptionOnFailure( ew < 1e-6 && eq < 1e-5 && eg < 1e-7 && es == 0.0,
                                 Marmot::Utilities::MakeString()
                                   << "tangent vs central differences at h = " << h << ": dw " << ew << ", dq/dh " << eq
                                   << ", dq/dgrad " << eg << ", dr/dh " << es );
      }
    }
    std::cout << "MoistureDiffusion: all tangents match central differences, max relative error " << maxErr
              << std::endl;
  }

  /// self-desiccation: the integral of the sink over consecutive increments is exact
  void testSelfDesiccation()
  {
    auto                                  material = create( fullBSB, 10 );
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    auto                                  fs = state( 0.9, 0.9, 0.0, 0.0, 0.0 );

    double consumed = 0.0, time = 0.0;
    for ( double dT : { 0.5, 3.0, 7.0, 20.0, 60.0 } ) {
      time += dT;
      evaluate( *material, fs, r, t, time, dT );
      throwExceptionOnFailure( r.source( 0 ) < 0.0, "self-desiccation is a sink" );
      consumed -= r.source( 0 ) * dT;
    }
    const double expected = 5.0 * ( 1.0 - std::exp( -time / 30.0 ) );
    throwExceptionOnFailure( checkIfEqual( consumed, expected, 1e-13 * expected ), "self-desiccation integral" );

    // increments much longer than the time constant (dT / tau_s > 709, where exp overflows): finite and still exact
    for ( double dT : { 25000.0, 1.0e7 } ) {
      time += dT;
      evaluate( *material, fs, r, t, time, dT );
      throwExceptionOnFailure( std::isfinite( r.source( 0 ) ) && r.source( 0 ) <= 0.0,
                               "self-desiccation sink of a very long increment" );
      consumed -= r.source( 0 ) * dT;
    }
    throwExceptionOnFailure( checkIfEqual( consumed, 5.0, 1e-13 * 5.0 ), "self-desiccation integral, long increments" );

    // one increment from the onset over (practically) the whole reaction
    evaluate( *material, fs, r, t, 1.0e6, 1.0e6 );
    throwExceptionOnFailure( checkIfEqual( -r.source( 0 ) * 1.0e6, 5.0, 1e-13 * 5.0 ),
                             "self-desiccation in one long increment" );
  }

  void testInvalidInput()
  {
    const auto throws = []( std::initializer_list< double > p ) {
      try {
        create( std::data( p ), static_cast< int >( p.size() ) );
      }
      catch ( const std::invalid_argument& ) {
        return true;
      }
      return false;
    };
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8 } ), "3 properties" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 0. } ), "5 properties" );
    throwExceptionOnFailure( throws( { 0., 0.05, 0.8, 15. } ), "D1 = 0" );
    throwExceptionOnFailure( throws( { 1., 0.0, 0.8, 15. } ), "alpha0 = 0" );
    throwExceptionOnFailure( throws( { 1., 1.5, 0.8, 15. } ), "alpha0 > 1" );
    throwExceptionOnFailure( throws( { 1., 0.05, 1.0, 15. } ), "hc = 1" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 0.5 } ), "n < 1" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 2., 1. } ), "isotherm type 2" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 0., -1. } ), "K < 0" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 1., 18. } ), "BSB with 6 properties" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 1., 18., 1.0, 20. } ), "BSB k = 1" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 0., 1., 0., 0., -1., 10. } ), "dwInf < 0" );
    throwExceptionOnFailure( throws( { 1., 0.05, 0.8, 15., 0., 1., 0., 0., 1., 0. } ), "tauS = 0" );
    throwExceptionOnFailure( !throws( { 1., 0.05, 0.8, 15., 0., 1., 0., 0., 0., 0. } ), "no sink, tauS unused" );
  }

} // namespace

int main()
{
  testClosedForms();
  testBSB();
  testTangentsByFiniteDifferences();
  testSelfDesiccation();
  testInvalidInput();
  return 0;
}
