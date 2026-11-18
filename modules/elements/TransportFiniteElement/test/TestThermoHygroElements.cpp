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
#include <vector>

/* The coupled heat and moisture transport elements DC*TH (fields "temperature", "relative humidity") with the
 * THERMOHYGRO material: consistent tangent of the element (all four blocks), of the latent heat film (temperature rows,
 * relative humidity columns), the uncoupled limit against the single field elements DC*T + DC*H (residual, tangent
 * and loads, to round-off), and zero time increments. */

using namespace Marmot::Testing;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Marmot::Discretizations::Elements::Element;

namespace {

  struct ElementType {
    std::string name; // without the field suffix
    int         nDim;
    int         nNodes;
    int         nFaces;
    bool        axisymmetric;
  };

  const std::vector< ElementType > elementTypes = {
    { "DC2D3", 2, 3, 3, false },
    { "DC2D4", 2, 4, 4, false },
    { "DC2D8", 2, 8, 4, false },
    { "DCAX4", 2, 4, 4, true },
    { "DCAX8", 2, 8, 4, true },
    { "DC3D4", 3, 4, 4, false },
    { "DC3D8", 3, 8, 6, false },
    { "DC3D20", 3, 20, 6, false },
  };

  const double thickness = 0.7;

  // heat (11), moisture (10, BSB isotherm and self-desiccation), couplings (11): all couplings on
  const std::vector< double > coupled = { 2.4,  0.9,  1.8,  1.6, 1.7,      0.1,  0.05, -0.08, 0.002, 0.001, 20.0,
                                          2e-3, 0.05, 0.8,  6.0, 1.0,      18.0, 0.73, 20.0,  5.0,   30.0,  -273.15,
                                          4700, 1.0,  0.01, 1.0, 4.186e-3, 0.3,  0.05, 4000,  60.0,  2.0 };
  // the same without couplings (21 properties)
  const std::vector< double > uncoupled( coupled.begin(), coupled.begin() + 21 );

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
    std::vector< double >      stateVars, initialStateVars;
    std::vector< double >      elementProperties{ thickness };
    int                        nDof;
  };

  TestElement makeElement( const std::string& name,
                           const ElementType& type,
                           const std::string& material,
                           const double*      p,
                           int                np,
                           bool               distort = true )
  {
    TestElement e;
    e.element.reset( Marmot::Registry::ElementFactory::createElement( name, 1 ) );
    e.coordinates = distortedCoordinates( type, distort );
    e.element->assignNodeCoordinates( e.coordinates.data() );
    e.element->assignProperty( Marmot::Discretizations::ElementProperties( e.elementProperties.data(), 1 ) );
    e.element->initializeYourself();
    e.element->assignProperty( Marmot::Discretizations::MaterialSection( material, p, np ) );
    e.stateVars.assign( e.element->getNumberOfRequiredStateVars(), 0.0 );
    e.element->assignStateVars( e.stateVars.data(), static_cast< int >( e.stateVars.size() ) );
    e.element->setInitialConditions( Element::MarmotMaterialInitialization, nullptr );
    // an equivalent age of 3 at the beginning of the increment
    for ( int q = 0; q < e.element->getNumberOfQuadraturePoints(); q++ )
      try {
        e.element->getStateView( "equivalent age", q ).stateLocation[0] = 3.0;
      }
      catch ( const std::exception& ) {
      }
    e.initialStateVars = e.stateVars;
    e.nDof             = e.element->getNDofPerElement();
    return e;
  }

  /// the kernels of an increment, always from the state at its beginning
  void kernels( TestElement&    e,
                const VectorXd& Q,
                const VectorXd& QOld,
                double          time,
                double          dT,
                VectorXd&       P,
                MatrixXd&       K )
  {
    e.stateVars       = e.initialStateVars;
    P                 = VectorXd::Zero( e.nDof );
    K                 = MatrixXd::Zero( e.nDof, e.nDof );
    const VectorXd dQ = Q - QOld;
    e.element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), time, dT );
  }

  /// nodal values: temperatures (dof 2 i) and relative humidities (dof 2 i + 1)
  VectorXd twoFields( int nNodes, double T0, double dT, double h0, double dh, double phase )
  {
    VectorXd Q( 2 * nNodes );
    for ( int i = 0; i < nNodes; i++ ) {
      Q( 2 * i )     = T0 + dT * std::sin( 1.3 * i + phase );
      Q( 2 * i + 1 ) = h0 + dh * std::cos( 0.9 * i + phase );
    }
    return Q;
  }

  /// element tangent against central differences, block by block
  void testTangent( const ElementType& type )
  {
    double maxErr = 0.0;
    for ( double dT : { 0.37, 0.0 } ) {
      auto     e    = makeElement( type.name + "TH", type, "THERMOHYGRO", coupled.data(), 32 );
      VectorXd Q    = twoFields( type.nNodes, 45.0, 15.0, 0.8, 0.12, 0.4 );
      VectorXd QOld = twoFields( type.nNodes, 30.0, 5.0, 0.9, 0.05, 1.1 );
      VectorXd P, Pp, Pm;
      MatrixXd K, Kd;
      kernels( e, Q, QOld, 5.0, dT, P, K );
      MatrixXd KFD( e.nDof, e.nDof );
      for ( int j = 0; j < e.nDof; j++ ) {
        const double h  = j % 2 == 0 ? 1e-4 : 1e-7;
        VectorXd     Qp = Q, Qm = Q;
        Qp( j ) += h;
        Qm( j ) -= h;
        kernels( e, Qp, QOld, 5.0, dT, Pp, Kd );
        kernels( e, Qm, QOld, 5.0, dT, Pm, Kd );
        KFD.col( j ) = ( Pp - Pm ) / ( 2 * h );
      }
      // the four blocks: rows / columns of T (even) and h (odd)
      for ( int a = 0; a < 2; a++ )
        for ( int b = 0; b < 2; b++ ) {
          MatrixXd Kab( type.nNodes, type.nNodes ), KFDab( type.nNodes, type.nNodes );
          for ( int i = 0; i < type.nNodes; i++ )
            for ( int j = 0; j < type.nNodes; j++ ) {
              Kab( i, j )   = K( 2 * i + a, 2 * j + b );
              KFDab( i, j ) = KFD( 2 * i + a, 2 * j + b );
            }
          const double err = ( Kab - KFDab ).norm() / KFDab.norm();
          maxErr           = std::max( maxErr, err );
          throwExceptionOnFailure( KFDab.norm() > 0.0 && err < 1e-6,
                                   Marmot::Utilities::MakeString()
                                     << type.name << "TH: tangent block (" << a << ", " << b << "), dT = " << dT
                                     << ": relative error " << err );
        }
    }
    std::cout << type.name
              << "TH: all four tangent blocks vs central differences (dT = 0.37 and 0), max relative error " << maxErr
              << std::endl;
  }

  /// latent heat film: temperature rows, relative humidity columns, on every face
  void testLatentHeatFilm( const ElementType& type )
  {
    auto           e      = makeElement( type.name + "TH", type, "THERMOHYGRO", coupled.data(), 32 );
    const VectorXd Q      = twoFields( type.nNodes, 45.0, 15.0, 0.8, 0.12, 0.4 );
    const double   film[] = { 25.0, 0.5 };
    for ( int face = 1; face <= type.nFaces; face++ ) {
      auto contribution = [&]( const VectorXd& Q_, VectorXd& R, MatrixXd& K ) {
        VectorXd Pext = VectorXd::Zero( e.nDof );
        K             = MatrixXd::Zero( e.nDof, e.nDof );
        e.element
          ->computeDistributedLoad( Element::LatentHeatFilm, Pext.data(), K.data(), face, film, Q_.data(), 1.0, 1.0 );
        R = -Pext;
      };
      VectorXd R, Rp, Rm;
      MatrixXd K, Kd, KFD( e.nDof, e.nDof );
      contribution( Q, R, K );
      for ( int j = 0; j < e.nDof; j++ ) {
        const double h  = 1e-6;
        VectorXd     Qp = Q, Qm = Q;
        Qp( j ) += h;
        Qm( j ) -= h;
        contribution( Qp, Rp, Kd );
        contribution( Qm, Rm, Kd );
        KFD.col( j ) = ( Rp - Rm ) / ( 2 * h );
      }
      bool onlyTRows = true, onlyHCols = true;
      for ( int i = 0; i < e.nDof; i++ ) {
        onlyTRows = onlyTRows && ( i % 2 == 0 || R( i ) == 0.0 );
        onlyHCols = onlyHCols && ( i % 2 == 1 || K.col( i ).norm() == 0.0 );
      }
      const double err = ( K - KFD ).norm() / KFD.norm();
      throwExceptionOnFailure( err < 1e-8 && onlyTRows && onlyHCols && R.sum() > 0.0,
                               type.name + "TH: latent heat film on face " + std::to_string( face ) );
    }
    // a single field element has no relative humidity to couple to
    auto     t = makeElement( type.name + "T", type, "HEATCONDUCTION", coupled.data(), 11 );
    VectorXd P = VectorXd::Zero( t.nDof ), QT = VectorXd::Constant( t.nDof, 20.0 );
    MatrixXd K     = MatrixXd::Zero( t.nDof, t.nDof );
    bool     threw = false;
    try {
      t.element->computeDistributedLoad( Element::LatentHeatFilm, P.data(), K.data(), 1, film, QT.data(), 1.0, 1.0 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, type.name + "T: a latent heat film needs the relative humidity" );
    std::cout << type.name << "TH: latent heat film tangent (T rows, h columns) on all faces" << std::endl;
  }

  /// uncoupled limit: DC*TH with THERMOHYGRO (21 properties) == DC*T with HEATCONDUCTION + DC*H with MOISTUREDIFFUSION
  void testUncoupledLimit( const ElementType& type )
  {
    auto th = makeElement( type.name + "TH", type, "THERMOHYGRO", uncoupled.data(), 21 );
    auto t  = makeElement( type.name + "T", type, "HEATCONDUCTION", coupled.data(), 11 );
    auto hy = makeElement( type.name + "H", type, "MOISTUREDIFFUSION", coupled.data() + 11, 10 );

    const VectorXd Q    = twoFields( type.nNodes, 45.0, 15.0, 0.8, 0.12, 0.4 ),
                   QOld = twoFields( type.nNodes, 30.0, 5.0, 0.9, 0.05, 1.1 );
    VectorXd QT( type.nNodes ), QTOld( type.nNodes ), QH( type.nNodes ), QHOld( type.nNodes );
    for ( int i = 0; i < type.nNodes; i++ ) {
      QT( i )    = Q( 2 * i );
      QTOld( i ) = QOld( 2 * i );
      QH( i )    = Q( 2 * i + 1 );
      QHOld( i ) = QOld( 2 * i + 1 );
    }
    double maxDiff = 0.0;
    auto   compare = [&]( const VectorXd&    P,
                        const MatrixXd&    K,
                        const VectorXd&    PT,
                        const MatrixXd&    KT,
                        const VectorXd&    PH,
                        const MatrixXd&    KH,
                        const std::string& what ) {
      double d = 0.0, n = 0.0;
      for ( int i = 0; i < type.nNodes; i++ ) {
        d = std::max( { d,
                        std::abs( P( 2 * i ) - PT( i ) ) / ( PT.cwiseAbs().maxCoeff() + 1e-300 ),
                        std::abs( P( 2 * i + 1 ) - PH( i ) ) / ( PH.cwiseAbs().maxCoeff() + 1e-300 ) } );
        for ( int j = 0; j < type.nNodes; j++ ) {
          d = std::max(
            { d,
              std::abs( K( 2 * i, 2 * j ) - KT( i, j ) ) / ( KT.cwiseAbs().maxCoeff() + 1e-300 ),
              std::abs( K( 2 * i + 1, 2 * j + 1 ) - KH( i, j ) ) / ( KH.cwiseAbs().maxCoeff() + 1e-300 ) } );
          n += std::abs( K( 2 * i, 2 * j + 1 ) ) + std::abs( K( 2 * i + 1, 2 * j ) );
        }
      }
      maxDiff = std::max( maxDiff, d );
      throwExceptionOnFailure( d < 1e-14 && n == 0.0,
                               Marmot::Utilities::MakeString() << type.name << "TH uncoupled, " << what
                                                               << ": difference " << d << ", coupling blocks " << n );
    };

    for ( double dT : { 0.37, 0.0 } ) {
      VectorXd P, PT, PH;
      MatrixXd K, KT, KH;
      kernels( th, Q, QOld, 5.0, dT, P, K );
      kernels( t, QT, QTOld, 5.0, dT, PT, KT );
      kernels( hy, QH, QHOld, 5.0, dT, PH, KH );
      compare( P, K, PT, KT, PH, KH, "kernels, dT = " + std::to_string( dT ) );
    }

    // loads: heat film and moisture film on face 1, body sources
    const double heatFilm[] = { 10.0, 20.0 }, moistureFilm[] = { 3e-2, 0.5 }, sources[] = { 5.0, -0.2 };
    VectorXd     P = VectorXd::Zero( th.nDof ), PT = VectorXd::Zero( t.nDof ), PH = VectorXd::Zero( hy.nDof );
    MatrixXd     K = MatrixXd::Zero( th.nDof, th.nDof ), KT = MatrixXd::Zero( t.nDof, t.nDof ),
             KH = MatrixXd::Zero( hy.nDof, hy.nDof );
    th.element->computeDistributedLoad( Element::HeatFilm, P.data(), K.data(), 1, heatFilm, Q.data(), 5.0, 0.37 );
    th.element
      ->computeDistributedLoad( Element::MoistureFilm, P.data(), K.data(), 1, moistureFilm, Q.data(), 5.0, 0.37 );
    th.element->computeBodyForce( P.data(), K.data(), sources, Q.data(), 5.0, 0.37 );
    t.element->computeDistributedLoad( Element::HeatFilm, PT.data(), KT.data(), 1, heatFilm, QT.data(), 5.0, 0.37 );
    t.element->computeBodyForce( PT.data(), KT.data(), &sources[0], QT.data(), 5.0, 0.37 );
    hy.element
      ->computeDistributedLoad( Element::MoistureFilm, PH.data(), KH.data(), 1, moistureFilm, QH.data(), 5.0, 0.37 );
    hy.element->computeBodyForce( PH.data(), KH.data(), &sources[1], QH.data(), 5.0, 0.37 );
    compare( P, K, PT, KT, PH, KH, "loads" );
    std::cout << type.name << "TH: uncoupled limit == " << type.name << "T + " << type.name
              << "H (kernels dT = 0.37 and 0, film and body loads), max relative difference " << maxDiff << std::endl;
  }

  /// zero time increment: no flow, no sources, no loads; the content increments balance
  void testZeroIncrement( const ElementType& type )
  {
    auto           e = makeElement( type.name + "TH", type, "THERMOHYGRO", coupled.data(), 32 );
    const VectorXd Q = twoFields( type.nNodes, 45.0, 15.0, 0.8, 0.12, 0.4 );
    VectorXd       P;
    MatrixXd       K;
    kernels( e, Q, Q, 5.0, 0.0, P, K );
    throwExceptionOnFailure( P.norm() == 0.0, type.name + "TH: nothing flows in a zero increment" );
    VectorXd     Pl     = VectorXd::Zero( e.nDof );
    const double film[] = { 25.0, 0.5 }, sources[] = { 1.0, 1.0 };
    e.element->computeDistributedLoad( Element::LatentHeatFilm, Pl.data(), K.data(), 1, film, Q.data(), 5.0, 0.0 );
    e.element->computeDistributedLoad( Element::HeatFilm, Pl.data(), K.data(), 1, film, Q.data(), 5.0, 0.0 );
    e.element->computeDistributedLoad( Element::MoistureFilm, Pl.data(), K.data(), 1, film, Q.data(), 5.0, 0.0 );
    e.element->computeBodyForce( Pl.data(), K.data(), sources, Q.data(), 5.0, 0.0 );
    throwExceptionOnFailure( Pl.norm() == 0.0, type.name + "TH: no loads in a zero increment" );

    // node fields
    const auto fields = e.element->getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == type.nNodes && fields[0].size() == 2 &&
                               fields[0][0] == "temperature" && fields[0][1] == "relative humidity",
                             type.name + "TH: fields temperature, relative humidity expected" );
    // a single field material is rejected
    bool threw = false;
    try {
      makeElement( type.name + "TH", type, "HEATCONDUCTION", coupled.data(), 11 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, type.name + "TH: a single field material must be rejected" );
  }

} // namespace

int main()
{
  for ( const auto& type : elementTypes ) {
    testTangent( type );
    testLatentHeatFilm( type );
    testUncoupledLimit( type );
    testZeroIncrement( type );
  }
  return 0;
}
