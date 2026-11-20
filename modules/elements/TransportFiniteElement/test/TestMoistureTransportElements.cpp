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
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/TransportFiniteElement.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

/* The moisture diffusion elements DC*H (field "relative humidity") with the MOISTUREDIFFUSION material: consistent
 * tangent of the element (Bazant-Najjar diffusivity, BSB isotherm, self-desiccation sink) and of the moisture film,
 * load resultants, zero time increments, steady state and the sealed (self-desiccating) balance. */

using namespace Marmot::Testing;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Marmot::Discretizations::Elements::Element;

namespace {

  struct ElementType {
    std::string name;
    int         nDim;
    int         nNodes;
    int         nFaces;
    bool        axisymmetric;
  };

  const std::vector< ElementType > elementTypes = {
    { "DC2D3H", 2, 3, 3, false },
    { "DC2D4H", 2, 4, 4, false },
    { "DC2D8H", 2, 8, 4, false },
    { "DCAX4H", 2, 4, 4, true },
    { "DCAX8H", 2, 8, 4, true },
    { "DC3D4H", 3, 4, 4, false },
    { "DC3D8H", 3, 8, 6, false },
    { "DC3D20H", 3, 20, 6, false },
  };

  const double thickness = 0.7;

  // D1, alpha0, hc, n, isotherm BSB (C, k, Vm), self-desiccation (dwInf, tauS)
  const double bsbMaterial[] = { 2.0e-3, 0.05, 0.8, 6.0, 1.0, 18.0, 0.73, 20.0, 5.0, 30.0 };
  // classical form: linear isotherm K = 1
  const double linearMaterial[] = { 2.0e-3, 0.05, 0.8, 6.0, 0.0, 1.0 };
  // steady state: K = 0
  const double steadyMaterial[] = { 2.0e-3, 0.05, 0.8, 6.0, 0.0, 0.0 };

  /// parametric node coordinates of the shapes, in the node order of Marmot
  std::vector< VectorXd > referenceNodes( int nDim, int nNodes )
  {
    std::vector< VectorXd > X;
    auto                    add = [&]( std::initializer_list< double > x ) {
      VectorXd v( nDim );
      int      i = 0;
      for ( double xi : x )
        v( i++ ) = xi;
      X.push_back( v );
    };
    if ( nDim == 2 && nNodes == 3 ) {
      add( { 0, 0 } );
      add( { 1, 0 } );
      add( { 0, 1 } );
    }
    else if ( nDim == 2 ) {
      add( { -1, -1 } );
      add( { 1, -1 } );
      add( { 1, 1 } );
      add( { -1, 1 } );
      if ( nNodes == 8 )
        for ( int i = 0; i < 4; i++ )
          X.push_back( 0.5 * ( X[i] + X[( i + 1 ) % 4] ) );
    }
    else if ( nNodes == 4 ) {
      add( { 0, 0, 0 } );
      add( { 1, 0, 0 } );
      add( { 0, 1, 0 } );
      add( { 0, 0, 1 } );
    }
    else {
      add( { -1, -1, -1 } );
      add( { 1, -1, -1 } );
      add( { 1, 1, -1 } );
      add( { -1, 1, -1 } );
      add( { -1, -1, 1 } );
      add( { 1, -1, 1 } );
      add( { 1, 1, 1 } );
      add( { -1, 1, 1 } );
      if ( nNodes == 20 ) {
        for ( int i = 0; i < 4; i++ )
          X.push_back( 0.5 * ( X[i] + X[( i + 1 ) % 4] ) );
        for ( int i = 0; i < 4; i++ )
          X.push_back( 0.5 * ( X[4 + i] + X[4 + ( i + 1 ) % 4] ) );
        for ( int i = 0; i < 4; i++ )
          X.push_back( 0.5 * ( X[i] + X[i + 4] ) );
      }
    }
    return X;
  }

  /// a multilinear distortion: straight edges stay straight and mid-side nodes stay at the middle of their edges
  VectorXd distortedCoordinates( const ElementType& type, bool distort = true )
  {
    const auto X = referenceNodes( type.nDim, type.nNodes );
    VectorXd   coords( type.nDim * type.nNodes );
    for ( int i = 0; i < type.nNodes; i++ ) {
      const VectorXd& x = X[i];
      VectorXd        y( type.nDim );
      if ( !distort )
        y = x;
      else if ( type.nDim == 2 ) {
        y( 0 ) = 1.2 * x( 0 ) + 0.3 * x( 1 ) + 0.15 * x( 0 ) * x( 1 );
        y( 1 ) = 0.1 * x( 0 ) + 0.9 * x( 1 ) - 0.1 * x( 0 ) * x( 1 );
      }
      else {
        y( 0 ) = 1.1 * x( 0 ) + 0.2 * x( 1 ) - 0.1 * x( 2 ) + 0.1 * x( 0 ) * x( 1 ) + 0.05 * x( 0 ) * x( 1 ) * x( 2 );
        y( 1 ) = 0.1 * x( 0 ) + 0.9 * x( 1 ) + 0.15 * x( 2 ) - 0.08 * x( 1 ) * x( 2 );
        y( 2 ) = -0.1 * x( 0 ) + 0.05 * x( 1 ) + 1.3 * x( 2 ) + 0.1 * x( 0 ) * x( 2 );
      }
      if ( type.axisymmetric )
        y( 0 ) += 2.5; // r > 0
      coords.segment( i * type.nDim, type.nDim ) = y;
    }
    return coords;
  }

  struct TestElement {
    std::unique_ptr< Element > element;
    VectorXd                   coordinates;
    std::vector< double >      stateVars;
    std::vector< double >      elementProperties{ thickness };
    int                        nDof;
  };

  TestElement makeElement( const ElementType& type, const double* p, int np, bool distort = true )
  {
    TestElement e;
    e.element.reset( Marmot::Registry::ElementFactory::createElement( type.name, 1 ) );
    e.coordinates = distortedCoordinates( type, distort );
    e.element->assignNodeCoordinates( e.coordinates.data() );
    e.element->assignProperty( Marmot::Discretizations::ElementProperties( e.elementProperties.data(), 1 ) );
    e.element->initializeYourself();
    e.element->assignProperty( Marmot::Discretizations::MaterialSection( "MOISTUREDIFFUSION", p, np ) );
    e.stateVars.assign( e.element->getNumberOfRequiredStateVars(), 0.0 );
    e.element->assignStateVars( e.stateVars.data(), static_cast< int >( e.stateVars.size() ) );
    e.nDof = e.element->getNDofPerElement();
    return e;
  }

  void kernels( TestElement&    e,
                const VectorXd& Q,
                const VectorXd& QOld,
                double          time,
                double          dT,
                VectorXd&       P,
                MatrixXd&       K )
  {
    P                 = VectorXd::Zero( e.nDof );
    K                 = MatrixXd::Zero( e.nDof, e.nDof );
    const VectorXd dQ = Q - QOld;
    e.element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), time, dT );
  }

  double measureBySource( TestElement& e )
  {
    VectorXd     P = VectorXd::Zero( e.nDof );
    MatrixXd     K = MatrixXd::Zero( e.nDof, e.nDof );
    const double s = 1.0;
    e.element->computeBodyForce( P.data(), K.data(), &s, P.data(), 1.0, 1.0 );
    return P.sum();
  }

  /// element tangent against central differences: nonlinear D(h), BSB isotherm, sink; h around h_c
  void testTangent( const ElementType& type )
  {
    double maxErr = 0.0;
    for ( auto [p, np] : { std::pair{ bsbMaterial, 10 }, std::pair{ linearMaterial, 6 } } ) {
      auto     e = makeElement( type, p, np );
      VectorXd Q( e.nDof ), QOld( e.nDof );
      for ( int i = 0; i < e.nDof; i++ ) {
        Q( i )    = 0.72 + 0.2 * std::sin( 1.3 * i + 0.4 );
        QOld( i ) = 0.85 + 0.1 * std::cos( 0.7 * i );
      }
      VectorXd P, Pp, Pm;
      MatrixXd K, Kd;
      kernels( e, Q, QOld, 5.0, 0.37, P, K );
      MatrixXd     KFD( e.nDof, e.nDof );
      const double h = 1e-6;
      for ( int j = 0; j < e.nDof; j++ ) {
        VectorXd Qp = Q, Qm = Q;
        Qp( j ) += h;
        Qm( j ) -= h;
        kernels( e, Qp, QOld, 5.0, 0.37, Pp, Kd );
        kernels( e, Qm, QOld, 5.0, 0.37, Pm, Kd );
        KFD.col( j ) = ( Pp - Pm ) / ( 2 * h );
      }
      const double err = ( K - KFD ).norm() / KFD.norm();
      maxErr           = std::max( maxErr, err );
      throwExceptionOnFailure( err < 1e-7, type.name + ": element tangent does not match finite differences" );
    }
    std::cout << type.name << ": element tangent vs central differences, relative error " << maxErr << std::endl;
  }

  /// moisture film tangent on every face, against central differences of the residual contribution -P
  void testFilmTangent( const ElementType& type )
  {
    auto     e = makeElement( type, bsbMaterial, 10 );
    VectorXd Q( e.nDof );
    for ( int i = 0; i < e.nDof; i++ )
      Q( i ) = 0.8 + 0.1 * std::sin( 1.1 * i );
    const double film[] = { 3.0e-2, 0.5 };
    for ( int face = 1; face <= type.nFaces; face++ ) {
      auto contribution = [&]( const VectorXd& Q_, VectorXd& R, MatrixXd& K ) {
        VectorXd Pext = VectorXd::Zero( e.nDof );
        K             = MatrixXd::Zero( e.nDof, e.nDof );
        e.element
          ->computeDistributedLoad( Element::MoistureFilm, Pext.data(), K.data(), face, film, Q_.data(), 1.0, 1.0 );
        R = -Pext;
      };
      VectorXd R, Rp, Rm;
      MatrixXd K, Kd, KFD( e.nDof, e.nDof );
      contribution( Q, R, K );
      const double h = 1e-6;
      for ( int j = 0; j < e.nDof; j++ ) {
        VectorXd Qp = Q, Qm = Q;
        Qp( j ) += h;
        Qm( j ) -= h;
        contribution( Qp, Rp, Kd );
        contribution( Qm, Rm, Kd );
        KFD.col( j ) = ( Rp - Rm ) / ( 2 * h );
      }
      const double err = ( K - KFD ).norm() / KFD.norm();
      throwExceptionOnFailure( err < 1e-8,
                               type.name + ": moisture film tangent on face " + std::to_string( face ) +
                                 " does not match finite differences" );
      throwExceptionOnFailure( R.sum() > 0.0, type.name + ": a wetter body must dry out" );
    }
    std::cout << type.name << ": moisture film tangent on all faces matches central differences" << std::endl;
  }

  /// moisture flux and film on face 1 of the undistorted element against flux times area; heat loads are rejected
  void testResultants( const ElementType& type )
  {
    auto   e    = makeElement( type, linearMaterial, 6, false );
    double area = 0.0;
    if ( type.nDim == 2 && type.nNodes == 3 )
      area = thickness;
    else if ( type.nDim == 3 && type.nNodes == 4 )
      area = 0.5;
    else if ( type.nDim == 3 )
      area = 4.0;
    else if ( type.axisymmetric )
      area = Marmot::Math::Constants::Pi * ( 3.5 * 3.5 - 1.5 * 1.5 );
    else
      area = 2.0 * thickness;

    VectorXd     P = VectorXd::Zero( e.nDof ), Q = VectorXd::Constant( e.nDof, 0.9 );
    MatrixXd     K = MatrixXd::Zero( e.nDof, e.nDof );
    const double q = 3.0e-3;
    e.element->computeDistributedLoad( Element::MoistureFlux, P.data(), K.data(), 1, &q, Q.data(), 1.0, 1.0 );
    throwExceptionOnFailure( checkIfEqual( P.sum(), q * area, 1e-12 * q * area ) && K.norm() == 0.0,
                             type.name + ": moisture flux resultant" );

    const double film[] = { 2.0e-2, 0.6 };
    P.setZero();
    e.element->computeDistributedLoad( Element::MoistureFilm, P.data(), K.data(), 1, film, Q.data(), 1.0, 1.0 );
    throwExceptionOnFailure( checkIfEqual( P.sum(), -2.0e-2 * 0.3 * area, 1e-12 * 6e-3 * area ),
                             type.name + ": moisture film resultant" );
    throwExceptionOnFailure( checkIfEqual( K.sum(), 2.0e-2 * area, 1e-12 * 2e-2 * area ),
                             type.name + ": moisture film tangent resultant" );

    // the element has no temperature: heat loads must be rejected
    for ( auto loadType : { Element::HeatFlux, Element::HeatFilm, Element::Pressure } ) {
      bool threw = false;
      try {
        e.element->computeDistributedLoad( loadType, P.data(), K.data(), 1, film, Q.data(), 1.0, 1.0 );
      }
      catch ( const std::invalid_argument& ) {
        threw = true;
      }
      throwExceptionOnFailure( threw, type.name + ": a heat or mechanical load must be rejected" );
    }
  }

  /// storage, the zero increment, the self-desiccation balance of a sealed element, and steady state
  void testStorageAndSink( const ElementType& type )
  {
    auto         e = makeElement( type, bsbMaterial, 10 );
    const double V = measureBySource( e );
    const auto   w = []( double h ) {
      return 18.0 * 0.73 * 20.0 * h / ( ( 1.0 - 0.73 * h ) * ( 1.0 + 17.0 * 0.73 * h ) );
    };

    // uniform drop 0.95 -> 0.9 over [time - dT, time] = [4, 6] in a sealed element: sum P = (dw / dT - r) V
    VectorXd     QOld = VectorXd::Constant( e.nDof, 0.95 ), Q = VectorXd::Constant( e.nDof, 0.9 );
    VectorXd     P;
    MatrixXd     K;
    const double time = 6.0, dT = 2.0;
    kernels( e, Q, QOld, time, dT, P, K );
    const double r        = -5.0 * ( std::exp( -4.0 / 30.0 ) - std::exp( -6.0 / 30.0 ) ) / dT;
    const double expected = ( ( w( 0.9 ) - w( 0.95 ) ) / dT - r ) * V;
    throwExceptionOnFailure( checkIfEqual( P.sum(), expected, 1e-12 * std::abs( expected ) ),
                             type.name + ": sealed balance with self-desiccation" );

    // zero time increment: P = content increment, no flux, no sink, no loads
    VectorXd QRamp( e.nDof );
    for ( int i = 0; i < e.nDof; i++ )
      QRamp( i ) = 0.6 + 0.3 * i / e.nDof;
    kernels( e, QRamp, QRamp, time, 0.0, P, K );
    throwExceptionOnFailure( P.norm() == 0.0, type.name + ": no flux and no sink in a zero time increment" );
    kernels( e, Q, QOld, time, 0.0, P, K );
    throwExceptionOnFailure( checkIfEqual( P.sum(), ( w( 0.9 ) - w( 0.95 ) ) * V, 1e-12 * V ),
                             type.name + ": zero time increment, content increment" );
    VectorXd     Pl = VectorXd::Zero( e.nDof );
    const double s = 1.0, film[] = { 2.0e-2, 0.5 }, q = 1.0;
    e.element->computeBodyForce( Pl.data(), K.data(), &s, Q.data(), time, 0.0 );
    e.element->computeDistributedLoad( Element::MoistureFilm, Pl.data(), K.data(), 1, film, Q.data(), time, 0.0 );
    e.element->computeDistributedLoad( Element::MoistureFlux, Pl.data(), K.data(), 1, &q, Q.data(), time, 0.0 );
    throwExceptionOnFailure( Pl.norm() == 0.0, type.name + ": no loads in a zero time increment" );

    // without storage the steady equations hold whatever the time increment, and the film acts also for dT = 0
    auto     steady = makeElement( type, steadyMaterial, 6 );
    VectorXd Ps, Ps0;
    MatrixXd Ks, Ks0;
    kernels( steady, QRamp, QOld, time, 0.0, Ps0, Ks0 );
    kernels( steady, QRamp, QOld, time, 0.5, Ps, Ks );
    throwExceptionOnFailure( Ps0.norm() > 0.0 && ( Ps - Ps0 ).norm() == 0.0 && ( Ks - Ks0 ).norm() == 0.0,
                             type.name + ": steady state must not depend on the time increment" );
    Pl.setZero();
    steady.element->computeDistributedLoad( Element::MoistureFilm, Pl.data(), K.data(), 1, film, Q.data(), time, 0.0 );
    throwExceptionOnFailure( Pl.norm() > 0.0, type.name + ": steady state: the film acts in a zero increment" );
    std::cout << type.name << ": storage, sealed self-desiccation balance, zero increment and steady state OK"
              << std::endl;
  }

} // namespace

int main()
{
  for ( const auto& type : elementTypes ) {
    const auto element = std::unique_ptr< Element >( Marmot::Registry::ElementFactory::createElement( type.name, 1 ) );
    const auto fields  = element->getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == type.nNodes && fields[0].size() == 1 &&
                               fields[0][0] == "relative humidity",
                             type.name + ": one field 'relative humidity' per node expected" );

    testTangent( type );
    testFilmTangent( type );
    testResultants( type );
    testStorageAndSink( type );
  }
  return 0;
}
