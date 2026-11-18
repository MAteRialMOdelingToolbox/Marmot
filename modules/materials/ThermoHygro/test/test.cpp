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
#include "Marmot/ThermoHygro.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

/* The coupled heat and moisture transport material THERMOHYGRO: uncoupled limit against HEATCONDUCTION +
 * MOISTUREDIFFUSION (exact), closed forms of the couplings, all tangent blocks by central differences, the
 * equivalent age, and rejected inputs. */

using namespace Marmot::Testing;
using Marmot::Materials::TransportMaterial;

namespace {

  std::unique_ptr< TransportMaterial > create( const std::string& name, const double* p, int n )
  {
    return std::unique_ptr< TransportMaterial >(
      Marmot::Registry::TransportMaterialFactory::createMaterial( name, p, n, 1 ) );
  }

  // heat: rho, c0, k11, k22, k33, k12, k13, k23, beta_k, beta_c, T_ref (degrees Celsius)
  // moisture: D1, alpha0, hc, n, BSB (C, k, Vm), dwInf, tauS
  // coupling: T_0, Q_D/R, m, kappa, chi, c_w, beta_kh, ell, E_h/R, Q_inf, tau_q
  const std::vector< double > coupled = { 2.4,  0.9,  1.8,  1.6, 1.7,      0.1,  0.05, -0.08, 0.002, 0.001, 20.0,
                                          2e-3, 0.05, 0.8,  6.0, 1.0,      18.0, 0.73, 20.0,  5.0,   30.0,  -273.15,
                                          4700, 1.0,  0.01, 1.0, 4.186e-3, 0.3,  0.05, 4000,  60.0,  2.0 };

  std::vector< double > uncoupledOf( const std::vector< double >& p )
  {
    std::vector< double > q( p.begin(), p.begin() + 21 );
    q.insert( q.end(), { -273.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } );
    return q;
  }

  struct Evaluation {
    TransportMaterial::Response           r;
    TransportMaterial::AlgorithmicTangent t;
    double                                te = 0.0;
  };

  /// evaluate with the equivalent age teN at the beginning of the increment
  Evaluation evaluate( const TransportMaterial&             material,
                       const TransportMaterial::FieldState& fs,
                       double                               time,
                       double                               dT,
                       double                               teN = 3.0 )
  {
    Evaluation            e;
    std::vector< double > sv( std::max( 1, material.getNumberOfRequiredStateVars() ), teN );
    material.initializeResponse( e.r, e.t );
    e.r.stateVars = sv.data();
    material.computeTransport( e.r, e.t, fs, { time, dT } );
    e.te = sv[0];
    return e;
  }

  TransportMaterial::FieldState twoFieldState( double T, double TOld, double h, double hOld )
  {
    TransportMaterial::FieldState fs;
    fs.field.resize( 2 );
    fs.fieldOld.resize( 2 );
    fs.field << T, h;
    fs.fieldOld << TOld, hOld;
    fs.gradient.resize( 6 );
    fs.gradient << 3.0, -1.5, 0.7, 0.4, -0.25, 0.6;
    return fs;
  }

  TransportMaterial::FieldState oneFieldState( double v, double vOld, const Eigen::Vector3d& grad )
  {
    TransportMaterial::FieldState fs;
    fs.field    = TransportMaterial::FieldVector::Constant( 1, v );
    fs.fieldOld = TransportMaterial::FieldVector::Constant( 1, vOld );
    fs.gradient = grad;
    return fs;
  }

  /// with all couplings off, both fields follow HEATCONDUCTION and MOISTUREDIFFUSION bit by bit
  void testUncoupledLimit()
  {
    const auto p21      = std::vector< double >( coupled.begin(), coupled.begin() + 21 );
    const auto p32      = uncoupledOf( coupled );
    auto       heat     = create( "HEATCONDUCTION", coupled.data(), 11 );
    auto       moisture = create( "MOISTUREDIFFUSION", coupled.data() + 11, 10 );
    for ( auto [p, np] : { std::pair{ p21.data(), 21 }, std::pair{ p32.data(), 32 } } ) {
      auto th = create( "THERMOHYGRO", p, np );
      throwExceptionOnFailure( th->getNumberOfFields() == 2 && th->hasStorage( 0 ) && th->hasStorage( 1 ),
                               "two fields with storage" );
      for ( double h : { 0.3, 0.75, 0.92, 1.0, 1.02 } )
        for ( double dT : { 0.0, 0.37 } ) {
          const auto fs = twoFieldState( 47.0, 31.0, h, 0.88 );
          const auto c  = evaluate( *th, fs, 5.0, dT );
          const auto eh = evaluate( *heat, oneFieldState( 47.0, 31.0, fs.gradient.head< 3 >() ), 5.0, dT );
          const auto em = evaluate( *moisture, oneFieldState( h, 0.88, fs.gradient.tail< 3 >() ), 5.0, dT );

          bool same = c.r.contentIncrement( 0 ) == eh.r.contentIncrement( 0 ) &&
                      c.r.contentIncrement( 1 ) == em.r.contentIncrement( 0 ) && c.r.flux.head< 3 >() == eh.r.flux &&
                      c.r.flux.tail< 3 >() == em.r.flux && c.r.source( 0 ) == eh.r.source( 0 ) &&
                      c.r.source( 1 ) == em.r.source( 0 );
          same = same && c.t.dContentIncrement_dField( 0, 0 ) == eh.t.dContentIncrement_dField( 0, 0 ) &&
                 c.t.dContentIncrement_dField( 1, 1 ) == em.t.dContentIncrement_dField( 0, 0 ) &&
                 c.t.dFlux_dField.block< 3, 1 >( 0, 0 ) == eh.t.dFlux_dField &&
                 c.t.dFlux_dField.block< 3, 1 >( 3, 1 ) == em.t.dFlux_dField &&
                 c.t.dFlux_dGradient.block< 3, 3 >( 0, 0 ) == eh.t.dFlux_dGradient &&
                 c.t.dFlux_dGradient.block< 3, 3 >( 3, 3 ) == em.t.dFlux_dGradient &&
                 c.t.dSource_dField( 1, 1 ) == em.t.dSource_dField( 0, 0 );
          // the coupling blocks vanish
          same = same && c.t.dContentIncrement_dField( 0, 1 ) == 0.0 && c.t.dContentIncrement_dField( 1, 0 ) == 0.0 &&
                 c.t.dFlux_dField.block< 3, 1 >( 0, 1 ).norm() == 0.0 &&
                 c.t.dFlux_dField.block< 3, 1 >( 3, 0 ).norm() == 0.0 &&
                 c.t.dFlux_dGradient.block< 3, 3 >( 0, 3 ).norm() == 0.0 &&
                 c.t.dFlux_dGradient.block< 3, 3 >( 3, 0 ).norm() == 0.0 && c.t.dSource_dField( 1, 0 ) == 0.0 &&
                 c.t.dSource_dField( 0, 0 ) == 0.0;
          throwExceptionOnFailure( same,
                                   Marmot::Utilities::MakeString()
                                     << np << " properties, h = " << h << ", dT = " << dT
                                     << ": the uncoupled limit must reproduce HEATCONDUCTION + MOISTUREDIFFUSION" );
        }
    }
    std::cout << "uncoupled limit (21 and 32 properties) == HEATCONDUCTION + MOISTUREDIFFUSION, bit by bit"
              << std::endl;
  }

  /// the couplings in closed form
  void testClosedForms()
  {
    auto         th = create( "THERMOHYGRO", coupled.data(), 32 );
    const double T = 60.0, TOld = 40.0, h = 0.85, hOld = 0.9;
    const auto   fs = twoFieldState( T, TOld, h, hOld );
    const auto   c  = evaluate( *th, fs, 5.0, 0.5, 3.0 );

    // hygrothermic coefficient: (dh/dT) at constant water content = kappa
    const double hygrothermic = -c.t.dContentIncrement_dField( 1, 0 ) / c.t.dContentIncrement_dField( 1, 1 );
    throwExceptionOnFailure( checkIfEqual( hygrothermic, 0.01, 1e-15 ), "hygrothermic coefficient" );

    // the moisture flux: D_BN(h) f(T) (grad h + chi h Theta_v / Ta^2 grad T)
    double                dD       = 0;
    auto                  moisture = create( "MOISTUREDIFFUSION", coupled.data() + 11, 10 );
    const auto            em  = evaluate( *moisture, oneFieldState( h, hOld, Eigen::Vector3d( 1, 0, 0 ) ), 5.0, 0.5 );
    const double          DBN = -em.r.flux( 0 );
    const double          Ta = T + 273.15, TaRef = 293.15;
    const double          f  = Ta / TaRef * std::exp( 4700.0 * ( 1.0 / TaRef - 1.0 / Ta ) );
    const double          g  = Marmot::Materials::ThermoHygro::vapourTemperature / ( Ta * Ta );
    const Eigen::Vector3d qW = -DBN * f * ( fs.gradient.tail< 3 >() + h * g * fs.gradient.head< 3 >() );
    throwExceptionOnFailure( ( c.r.flux.tail< 3 >() - qW ).norm() < 1e-13 * qW.norm(), "moisture flux" );
    (void)dD;

    // the heat flux: (1 + beta_kh h) Fourier + ell q_w
    auto                  heat = create( "HEATCONDUCTION", coupled.data(), 11 );
    const auto            eh   = evaluate( *heat, oneFieldState( T, TOld, fs.gradient.head< 3 >() ), 5.0, 0.5 );
    const Eigen::Vector3d qT   = ( 1.0 + 0.3 * h ) * eh.r.flux + 0.05 * qW;
    throwExceptionOnFailure( ( c.r.flux.head< 3 >() - qT ).norm() < 1e-13 * qT.norm(), "heat flux" );

    // enthalpy increment: dry + c_w w (T - T_n)
    const double w   = c.r.contentIncrement( 1 ) + 0.0; // increment; recompute w from the isotherm instead
    const auto   bsb = []( double x ) {
      return 18.0 * 0.73 * 20.0 * x / ( ( 1.0 - 0.73 * x ) * ( 1.0 + 17.0 * 0.73 * x ) );
    };
    (void)w;
    const double wEnd = bsb( h - 0.01 * ( T - 20.0 ) ), wOld = bsb( hOld - 0.01 * ( TOld - 20.0 ) );
    throwExceptionOnFailure( checkIfEqual( c.r.contentIncrement( 1 ), wEnd - wOld, 1e-12 * wEnd ), "water content" );
    const double dH = eh.r.contentIncrement( 0 ) + 4.186e-3 * wEnd * ( T - TOld );
    throwExceptionOnFailure( checkIfEqual( c.r.contentIncrement( 0 ), dH, 1e-13 * std::abs( dH ) ), "enthalpy" );

    // equivalent age and sources: t_e = t_e,n + a(T) dT, exact averages of the exponential laws
    const double a  = std::exp( 4000.0 * ( 1.0 / TaRef - 1.0 / Ta ) );
    const double te = 3.0 + a * 0.5;
    throwExceptionOnFailure( checkIfEqual( c.te, te, 1e-14 * te ), "equivalent age" );
    const double rT = 60.0 * ( std::exp( -3.0 / 2.0 ) - std::exp( -te / 2.0 ) ) / 0.5;
    const double rH = -5.0 * ( std::exp( -3.0 / 30.0 ) - std::exp( -te / 30.0 ) ) / 0.5;
    throwExceptionOnFailure( checkIfEqual( c.r.source( 0 ), rT, 1e-12 * rT ), "heat of hydration" );
    throwExceptionOnFailure( checkIfEqual( c.r.source( 1 ), rH, -1e-12 * rH ), "self-desiccation" );
    std::cout << "closed forms of the couplings OK (maturity rate at 60 C: " << a << ")" << std::endl;
  }

  /// all four blocks of all tangents by central differences, for several states and both time increments
  void testTangentsByFiniteDifferences()
  {
    auto   th     = create( "THERMOHYGRO", coupled.data(), 32 );
    double maxErr = 0.0;
    for ( double T : { 5.0, 20.0, 65.0 } )
      for ( double h : { 0.4, 0.78, 0.83, 0.97 } )
        for ( double dT : { 0.0, 0.37 } ) {
          const auto fs = twoFieldState( T, T - 12.0, h, 0.9 );
          const auto c  = evaluate( *th, fs, 5.0, dT );

          const double scale[2] = { 1e-4, 1e-7 }; // perturbation of T, h
          for ( int b = 0; b < 2; b++ ) {
            auto fp = fs, fm = fs;
            fp.field( b ) += scale[b];
            fm.field( b ) -= scale[b];
            const auto ep = evaluate( *th, fp, 5.0, dT ), em = evaluate( *th, fm, 5.0, dT );
            for ( int a = 0; a < 2; a++ ) {
              const double dm  = ( ep.r.contentIncrement( a ) - em.r.contentIncrement( a ) ) / ( 2 * scale[b] );
              const double ds  = ( ep.r.source( a ) - em.r.source( a ) ) / ( 2 * scale[b] );
              const double nm  = c.t.dContentIncrement_dField.row( a ).cwiseAbs().maxCoeff() + 1e-30;
              const double ns  = c.t.dSource_dField.row( a ).cwiseAbs().maxCoeff() + 1e-30;
              const double em_ = std::abs( dm - c.t.dContentIncrement_dField( a, b ) ) / nm;
              const double es_ = ds == 0.0 && c.t.dSource_dField( a, b ) == 0.0
                                   ? 0.0
                                   : std::abs( ds - c.t.dSource_dField( a, b ) ) / ns;
              double       eq_ = 0.0;
              const double nq  = c.t.dFlux_dField.block( 3 * a, 0, 3, 2 ).cwiseAbs().maxCoeff() + 1e-30;
              for ( int i = 0; i < 3; i++ ) {
                const double dq = ( ep.r.flux( 3 * a + i ) - em.r.flux( 3 * a + i ) ) / ( 2 * scale[b] );
                eq_             = std::max( eq_, std::abs( dq - c.t.dFlux_dField( 3 * a + i, b ) ) / nq );
              }
              maxErr = std::max( { maxErr, em_, es_, eq_ } );
              throwExceptionOnFailure( em_ < 1e-6 && es_ < 1e-6 && eq_ < 1e-6,
                                       Marmot::Utilities::MakeString()
                                         << "field tangent block (" << a << ", " << b << ") at T = " << T
                                         << ", h = " << h << ", dT = " << dT << ": content " << em_ << ", source "
                                         << es_ << ", flux " << eq_ );
            }
          }
          for ( int j = 0; j < 6; j++ ) {
            const double e  = 1e-5;
            auto         fp = fs, fm = fs;
            fp.gradient( j ) += e;
            fm.gradient( j ) -= e;
            const auto ep = evaluate( *th, fp, 5.0, dT ), em = evaluate( *th, fm, 5.0, dT );
            for ( int i = 0; i < 6; i++ ) {
              const double dq = ( ep.r.flux( i ) - em.r.flux( i ) ) / ( 2 * e );
              const double n  = c.t.dFlux_dGradient.row( i ).cwiseAbs().maxCoeff() + 1e-30;
              const double er = std::abs( dq - c.t.dFlux_dGradient( i, j ) ) / n;
              maxErr          = std::max( maxErr, er );
              throwExceptionOnFailure( er < 1e-7, "gradient tangent" );
            }
          }
        }
    std::cout << "all tangent blocks vs central differences, max relative error " << maxErr << std::endl;
  }

  /// storage flags and rejected inputs
  void testStorageAndInputs()
  {
    auto p = coupled;
    p[0]   = 0.0; // rho = 0, but water with heat capacity: storage of heat
    throwExceptionOnFailure( create( "THERMOHYGRO", p.data(), 32 )->hasStorage( 0 ), "c_w w: heat storage" );
    p[26] = 0.0;
    throwExceptionOnFailure( !create( "THERMOHYGRO", p.data(), 32 )->hasStorage( 0 ), "no heat storage" );
    auto q = coupled;
    q[15]  = 0.0;
    q[16]  = 0.0; // linear isotherm, K = 0
    auto m = create( "THERMOHYGRO", q.data(), 32 );
    throwExceptionOnFailure( m->hasStorage( 0 ) && !m->hasStorage( 1 ), "steady moisture" );

    int rejected = 0;
    for ( auto [index, value] : std::vector< std::pair< int, double > >{ { 22, -1.0 },
                                                                         { 25, -0.5 },
                                                                         { 26, -1.0 },
                                                                         { 27, -1.0 },
                                                                         { 29, -1.0 },
                                                                         { 30, -1.0 },
                                                                         { 31, 0.0 },
                                                                         { 21, 30.0 },
                                                                         { 11, 0.0 },
                                                                         { 0, -1.0 } } ) {
      auto bad   = coupled;
      bad[index] = value;
      try {
        create( "THERMOHYGRO", bad.data(), 32 );
      }
      catch ( const std::invalid_argument& ) {
        rejected++;
      }
    }
    for ( int n : { 11, 20, 22, 31, 33 } )
      try {
        create( "THERMOHYGRO", coupled.data(), n );
      }
      catch ( const std::invalid_argument& ) {
        rejected++;
      }
    throwExceptionOnFailure( rejected == 15, Marmot::Utilities::MakeString() << "rejected " << rejected << " of 15" );

    // below absolute zero: reported
    bool       threw = false;
    const auto th    = create( "THERMOHYGRO", coupled.data(), 32 );
    try {
      evaluate( *th, twoFieldState( -300.0, 20.0, 0.8, 0.8 ), 5.0, 0.5 );
    }
    catch ( const std::runtime_error& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, "absolute temperature <= 0 must throw" );
    std::cout << "storage flags and 15 invalid inputs OK" << std::endl;
  }

} // namespace

int main()
{
  testUncoupledLimit();
  testClosedForms();
  testTangentsByFiniteDifferences();
  testStorageAndInputs();
  return 0;
}
