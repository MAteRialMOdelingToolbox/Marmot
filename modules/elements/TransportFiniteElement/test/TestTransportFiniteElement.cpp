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
#include "Marmot/MarmotTransportMaterial.h"
#include "Marmot/MarmotTransportMaterialFactory.h"
#include "Marmot/TransportFiniteElement.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

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
    { "DC2D3T", 2, 3, 3, false },
    { "DC2D4T", 2, 4, 4, false },
    { "DC2D8T", 2, 8, 4, false },
    { "DCAX4T", 2, 4, 4, true },
    { "DCAX8T", 2, 8, 4, true },
    { "DC3D4T", 3, 4, 4, false },
    { "DC3D8T", 3, 8, 6, false },
    { "DC3D20T", 3, 20, 6, false },
  };

  const double thickness = 0.7;

  // the anisotropic, temperature dependent material: rho, c, k11, k22, k33, k12, k13, k23, beta_k, beta_c, T_ref
  const double nonlinearMaterial[] = { 2.0, 3.0, 2.5, 1.5, 1.0, 0.3, -0.2, 0.1, 0.004, 0.002, 20.0 };
  // isotropic constant: rho, c, k
  const double linearIsotropicMaterial[] = { 2.0, 3.0, 1.7 };
  // anisotropic constant
  const double linearAnisotropicMaterial[] = { 2.0, 3.0, 2.5, 1.5, 1.0, 0.3, -0.2, 0.1, 0.0, 0.0, 0.0 };

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

  TestElement makeElement( const ElementType& type,
                           const double*      materialProperties,
                           int                nMaterialProperties,
                           bool               distort = true )
  {
    TestElement e;
    e.element.reset( Marmot::Registry::ElementFactory::createElement( type.name, 1 ) );
    e.coordinates = distortedCoordinates( type, distort );
    e.element->assignNodeCoordinates( e.coordinates.data() );
    e.element->assignProperty( Marmot::Discretizations::ElementProperties( e.elementProperties.data(), 1 ) );
    e.element->initializeYourself();
    e.element->assignProperty(
      Marmot::Discretizations::MaterialSection( "HEATCONDUCTION", materialProperties, nMaterialProperties ) );
    e.stateVars.assign( e.element->getNumberOfRequiredStateVars(), 0.0 );
    e.element->assignStateVars( e.stateVars.data(), static_cast< int >( e.stateVars.size() ) );
    e.nDof = e.element->getNDofPerElement();
    return e;
  }

  /// residual and tangent for the end-of-increment fields Q at fixed start-of-increment fields QOld
  void kernels( TestElement& e, const VectorXd& Q, const VectorXd& QOld, double dT, VectorXd& P, MatrixXd& K )
  {
    P                 = VectorXd::Zero( e.nDof );
    K                 = MatrixXd::Zero( e.nDof, e.nDof );
    const VectorXd dQ = Q - QOld;
    e.element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0.0, dT );
  }

  double relativeError( const MatrixXd& A, const MatrixXd& B )
  {
    return ( A - B ).norm() / B.norm();
  }

  /// the tangent of computeKernels against central differences, with a temperature dependent anisotropic material
  void testTangent( const ElementType& type )
  {
    auto e = makeElement( type, nonlinearMaterial, 11 );

    VectorXd Q( e.nDof ), QOld( e.nDof );
    for ( int i = 0; i < e.nDof; i++ ) {
      Q( i )    = 40.0 + 15.0 * std::sin( 1.3 * i + 0.4 );
      QOld( i ) = 30.0 + 10.0 * std::cos( 0.7 * i );
    }

    VectorXd P, Pp, Pm;
    MatrixXd K, Kd;
    kernels( e, Q, QOld, 0.37, P, K );

    MatrixXd     KFD( e.nDof, e.nDof );
    const double h = 1e-5;
    for ( int j = 0; j < e.nDof; j++ ) {
      VectorXd Qp = Q, Qm = Q;
      Qp( j ) += h;
      Qm( j ) -= h;
      kernels( e, Qp, QOld, 0.37, Pp, Kd );
      kernels( e, Qm, QOld, 0.37, Pm, Kd );
      KFD.col( j ) = ( Pp - Pm ) / ( 2 * h );
    }
    const double err = relativeError( K, KFD );
    std::cout << type.name << ": element tangent vs central differences, relative error " << err << std::endl;
    throwExceptionOnFailure( err < 1e-8, type.name + ": element tangent does not match finite differences" );

    // the state written to the quadrature points
    const auto flux = e.element->getStateView( "flux", 0 );
    throwExceptionOnFailure( flux.stateSize == 3, type.name + ": flux state of size 3 expected" );
  }

  /// the film and radiation tangents against central differences of the residual contribution -P
  void testLoadTangents( const ElementType& type )
  {
    auto e = makeElement( type, nonlinearMaterial, 11 );

    VectorXd Q( e.nDof );
    for ( int i = 0; i < e.nDof; i++ )
      Q( i ) = 400.0 + 50.0 * std::sin( 1.1 * i );

    const double film[]      = { 7.0, 290.0 };
    const double radiation[] = { 0.8 * 5.67e-8, 300.0, -273.15 };

    for ( auto [loadType, load, label] : { std::tuple{ Element::HeatFilm, film, "film" },
                                           std::tuple{ Element::HeatRadiation, radiation, "radiation" } } )
      for ( int face = 1; face <= type.nFaces; face++ ) {
        auto contribution = [&]( const VectorXd& Q_, VectorXd& R, MatrixXd& K ) {
          VectorXd Pext = VectorXd::Zero( e.nDof );
          K             = MatrixXd::Zero( e.nDof, e.nDof );
          e.element->computeDistributedLoad( loadType, Pext.data(), K.data(), face, load, Q_.data(), 0.0, 1.0 );
          R = -Pext;
        };
        VectorXd R, Rp, Rm;
        MatrixXd K, Kd, KFD( e.nDof, e.nDof );
        contribution( Q, R, K );
        const double h = 1e-4;
        for ( int j = 0; j < e.nDof; j++ ) {
          VectorXd Qp = Q, Qm = Q;
          Qp( j ) += h;
          Qm( j ) -= h;
          contribution( Qp, Rp, Kd );
          contribution( Qm, Rm, Kd );
          KFD.col( j ) = ( Rp - Rm ) / ( 2 * h );
        }
        const double err = relativeError( K, KFD );
        throwExceptionOnFailure( err < 1e-7,
                                 type.name + ": " + label + " tangent on face " + std::to_string( face ) +
                                   " does not match finite differences, relative error " + std::to_string( err ) );
        // heat flows out of a body warmer than its surroundings
        throwExceptionOnFailure( R.sum() > 0.0, type.name + ": " + label + " must remove heat" );
      }
    std::cout << type.name << ": film and radiation tangents on all faces match central differences" << std::endl;
  }

  /// the measure of the element (area times thickness, ring volume, or volume) by the body source
  double measureBySource( TestElement& e )
  {
    VectorXd     P = VectorXd::Zero( e.nDof );
    MatrixXd     K = MatrixXd::Zero( e.nDof, e.nDof );
    const double s = 1.0;
    e.element->computeBodyForce( P.data(), K.data(), &s, P.data(), 0.0, 1.0 );
    return P.sum();
  }

  /// flux and source loads on the undistorted reference element, against their resultants
  void testResultants( const ElementType& type )
  {
    auto e = makeElement( type, linearIsotropicMaterial, 3, false );

    // body source: the element measure
    double measure = 0.0;
    if ( type.nDim == 2 && type.nNodes == 3 )
      measure = 0.5 * thickness;
    else if ( type.nDim == 3 && type.nNodes == 4 )
      measure = 1.0 / 6.0;
    else if ( type.nDim == 3 )
      measure = 8.0;
    else if ( type.axisymmetric )
      measure = Marmot::Math::Constants::Pi * ( 3.5 * 3.5 - 1.5 * 1.5 ) * 2.0; // r in [1.5, 3.5], z in [-1, 1]
    else
      measure = 4.0 * thickness;
    throwExceptionOnFailure( checkIfEqual( measureBySource( e ), measure, 1e-12 * measure ),
                             type.name + ": body source resultant" );

    // surface flux on face 1: its area (times thickness, or over the ring)
    double area = 0.0;
    if ( type.nDim == 2 && type.nNodes == 3 )
      area = thickness;
    else if ( type.nDim == 3 && type.nNodes == 4 )
      area = 0.5;
    else if ( type.nDim == 3 )
      area = 4.0;
    else if ( type.axisymmetric )
      area = Marmot::Math::Constants::Pi * ( 3.5 * 3.5 - 1.5 * 1.5 ); // face 1 is the bottom, z = -1
    else
      area = 2.0 * thickness;
    VectorXd     P = VectorXd::Zero( e.nDof ), Q = VectorXd::Zero( e.nDof );
    MatrixXd     K = MatrixXd::Zero( e.nDof, e.nDof );
    const double q = 3.0;
    e.element->computeDistributedLoad( Element::HeatFlux, P.data(), K.data(), 1, &q, Q.data(), 0.0, 1.0 );
    throwExceptionOnFailure( checkIfEqual( P.sum(), q * area, 1e-12 * q * area ), type.name + ": heat flux resultant" );
    throwExceptionOnFailure( K.norm() == 0.0, type.name + ": a prescribed heat flux has no tangent" );

    // a uniform film flux h (T - T_inf) equals h (T - T_inf) times the area
    const double film[] = { 2.0, 10.0 };
    Q.setConstant( 25.0 );
    P.setZero();
    e.element->computeDistributedLoad( Element::HeatFilm, P.data(), K.data(), 1, film, Q.data(), 0.0, 1.0 );
    throwExceptionOnFailure( checkIfEqual( P.sum(), -2.0 * 15.0 * area, 1e-12 * 30.0 * area ),
                             type.name + ": film resultant" );

    // a mechanical load type is not a transport load
    bool threw = false;
    try {
      e.element->computeDistributedLoad( Element::Pressure, P.data(), K.data(), 1, &q, Q.data(), 0.0, 1.0 );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw, type.name + ": a pressure must be rejected" );
    std::cout << type.name << ": body source, heat flux and film resultants are exact" << std::endl;
  }

  /// steady state, linear temperature field, constant conductivity: the element residual is the negative outward
  /// heat flow through its faces, -sum_faces int N q.n dGamma, on a distorted element (discrete divergence theorem)
  void testLinearFieldConsistency( const ElementType& type )
  {
    // axisymmetric: a gradient along z with isotropic k only, else the flux is not divergence free in the ring
    auto e = type.axisymmetric ? makeElement( type, linearIsotropicMaterial, 3 )
                               : makeElement( type, linearAnisotropicMaterial, 11 );

    Eigen::Vector3d gradT( 0.7, -1.3, 0.4 );
    if ( type.nDim == 2 )
      gradT( 2 ) = 0.0;
    if ( type.axisymmetric )
      gradT << 0.0, -1.3, 0.0;

    Eigen::Matrix3d k;
    if ( type.axisymmetric )
      k = 1.7 * Eigen::Matrix3d::Identity();
    else
      k << 2.5, 0.3, -0.2, 0.3, 1.5, 0.1, -0.2, 0.1, 1.0;
    const Eigen::Vector3d q3 = -k * gradT;
    // dynamic copies of the in-plane parts: GCC reports false array-bounds warnings on head( n ) of a fixed size vector
    const VectorXd q     = q3.head( type.nDim );
    const VectorXd gradX = gradT.head( type.nDim );

    VectorXd Q( e.nDof );
    for ( int i = 0; i < type.nNodes; i++ )
      Q( i ) = 5.0 + gradX.dot( e.coordinates.segment( i * type.nDim, type.nDim ) );

    VectorXd P;
    MatrixXd K;
    kernels( e, Q, Q, 1.0, P, K );

    VectorXd outflow = VectorXd::Zero( e.nDof );
    for ( int face = 1; face <= type.nFaces; face++ ) {
      const auto shape = Marmot::Discretizations::Elements::getElementShapeByMetric( type.nDim, type.nNodes );
      Marmot::Discretizations::Elements::BoundaryElement be( shape, face, type.nDim, e.coordinates, type.axisymmetric );
      const VectorXd         Nn    = be.computeSurfaceNormalVectorialLoadVector(); // int N_i n dGamma
      const Eigen::VectorXi& nodes = be.getBoundaryNodeIndicesInParent();
      for ( int k_ = 0; k_ < nodes.size(); k_++ )
        outflow( nodes( k_ ) ) += Nn.segment( k_ * type.nDim, type.nDim ).dot( q );
    }
    if ( type.nDim == 2 && !type.axisymmetric )
      outflow *= thickness;

    const double err = ( P + outflow ).norm() / outflow.norm();
    std::cout << type.name << ": linear field, residual vs boundary heat flow, relative error " << err << std::endl;
    throwExceptionOnFailure( err < 1e-12, type.name + ": linear field consistency failed" );
    throwExceptionOnFailure( std::abs( P.sum() ) < 1e-12 * P.norm(), type.name + ": steady balance not conserved" );
  }

  /// transient term: a uniform temperature rise dT_ of the whole element stores rho c V dT_ / dt
  void testStorage( const ElementType& type )
  {
    auto     e    = makeElement( type, linearIsotropicMaterial, 3 );
    VectorXd QOld = VectorXd::Constant( e.nDof, 20.0 ), Q = VectorXd::Constant( e.nDof, 23.0 );
    VectorXd P;
    MatrixXd K;
    kernels( e, Q, QOld, 0.5, P, K );
    const double expected = 2.0 * 3.0 * measureBySource( e ) * 3.0 / 0.5;
    throwExceptionOnFailure( checkIfEqual( P.sum(), expected, 1e-12 * expected ), type.name + ": storage term" );
    // and the capacity matrix sums to rho c V / dt
    throwExceptionOnFailure( checkIfEqual( K.sum(), expected / 3.0, 1e-12 * expected ), type.name + ": capacity" );

    // a zero time increment lets no time pass: the content increment itself is balanced, nothing flows, and no load
    // acts; a nonuniform field shows that the flux is off
    VectorXd QRamp( e.nDof );
    for ( int i = 0; i < e.nDof; i++ )
      QRamp( i ) = 23.0 + 2.0 * i;
    kernels( e, QRamp, QRamp, 0.0, P, K );
    throwExceptionOnFailure( P.norm() == 0.0, type.name + ": no flux in a zero time increment" );
    kernels( e, Q, QOld, 0.0, P, K );
    throwExceptionOnFailure( checkIfEqual( P.sum(), expected * 0.5, 1e-12 * expected ),
                             type.name + ": zero time increment, content increment" );
    throwExceptionOnFailure( checkIfEqual( K.sum(), expected * 0.5 / 3.0, 1e-12 * expected ),
                             type.name + ": zero time increment, capacity" );
    VectorXd     Pl = VectorXd::Zero( e.nDof );
    const double s = 1.0, film[] = { 2.0, 10.0 };
    e.element->computeBodyForce( Pl.data(), K.data(), &s, Q.data(), 0.0, 0.0 );
    e.element->computeDistributedLoad( Element::HeatFilm, Pl.data(), K.data(), 1, film, Q.data(), 0.0, 0.0 );
    throwExceptionOnFailure( Pl.norm() == 0.0, type.name + ": no loads in a zero time increment" );

    // without storage (zero heat capacity) the steady state equations hold also for a zero time increment
    const double steadyMaterial[] = { 0.0, 3.0, 1.7 };
    auto         steady           = makeElement( type, steadyMaterial, 3 );
    VectorXd     Ps, Ps0;
    MatrixXd     Ks, Ks0;
    kernels( steady, QRamp, QOld, 0.0, Ps0, Ks0 );
    kernels( steady, QRamp, QOld, 0.5, Ps, Ks );
    throwExceptionOnFailure( Ps0.norm() > 0.0 && ( Ps - Ps0 ).norm() == 0.0 && ( Ks - Ks0 ).norm() == 0.0,
                             type.name + ": steady state must not depend on the time increment" );
  }

  /// whether a call throws (any std::exception)
  bool throws( const std::function< void() >& f )
  {
    try {
      f();
    }
    catch ( const std::exception& ) {
      return true;
    }
    return false;
  }

  /// the interface the solver uses besides the kernels: dof layout, shape, coordinates, initialization, states
  void testInterface( const ElementType& type )
  {
    auto e = makeElement( type, linearIsotropicMaterial, 3, false );

    const auto pattern = e.element->getDofIndicesPermutationPattern();
    throwExceptionOnFailure( static_cast< int >( pattern.size() ) == e.nDof, type.name + ": permutation size" );
    for ( int i = 0; i < e.nDof; i++ )
      throwExceptionOnFailure( pattern[i] == i, type.name + ": the dofs are ordered node by node already" );
    throwExceptionOnFailure( e.element->getNNodes() == type.nNodes && e.element->getNSpatialDimensions() == type.nDim,
                             type.name + ": number of nodes or dimensions" );
    throwExceptionOnFailure( !e.element->getElementShape().empty(), type.name + ": element shape" );

    // on the undistorted element the center is the mean of the nodes
    const auto center = e.element->getCoordinatesAtCenter();
    for ( int d = 0; d < type.nDim; d++ ) {
      double mean = 0.0;
      for ( int i = 0; i < type.nNodes; i++ )
        mean += e.coordinates( i * type.nDim + d ) / type.nNodes;
      throwExceptionOnFailure( checkIfEqual( center[d], mean, 1e-12 ), type.name + ": element center" );
    }

    // the quadrature points lie inside the element's bounding box
    const auto qpCoords = e.element->getCoordinatesAtQuadraturePoints();
    throwExceptionOnFailure( static_cast< int >( qpCoords.size() ) == e.element->getNumberOfQuadraturePoints(),
                             type.name + ": one coordinate per quadrature point" );
    for ( const auto& x : qpCoords )
      for ( int d = 0; d < type.nDim; d++ ) {
        double lo = 1e300, hi = -1e300;
        for ( int i = 0; i < type.nNodes; i++ ) {
          lo = std::min( lo, e.coordinates( i * type.nDim + d ) );
          hi = std::max( hi, e.coordinates( i * type.nDim + d ) );
        }
        throwExceptionOnFailure( x[d] > lo && x[d] < hi, type.name + ": quadrature point outside the element" );
      }

    // material initialization is the only initial condition; the element's own states and the material's are viewable
    e.element->setInitialConditions( Element::MarmotMaterialInitialization, nullptr );
    throwExceptionOnFailure( throws( [&] { e.element->setInitialConditions( Element::Sigma11, nullptr ); } ),
                             type.name + ": a stress initial condition must be rejected" );
    throwExceptionOnFailure( e.element->getStateView( "flux", 0 ).stateSize == 3, type.name + ": flux state" );
    throwExceptionOnFailure( throws( [&] { e.element->getStateView( "no such state", 0 ); } ),
                             type.name + ": an unknown state must be rejected" );
  }

  /// a transport material that only states its number of fields: the defaults of the base class
  class MinimalTransportMaterial : public Marmot::Materials::TransportMaterial {
  public:
    using TransportMaterial::TransportMaterial;
    int  getNumberOfFields() const override { return 1; }
    void computeTransport( Response&, AlgorithmicTangent&, const FieldState&, const TimeIncrement& ) const override {}
  };

  /// invalid input is rejected with an exception
  void testInvalidInput()
  {
    using Marmot::Discretizations::Elements::TransportFiniteElement;
    namespace Quadrature = Marmot::Discretizations::Elements::Quadrature;

    throwExceptionOnFailure( throws( [] {
                               Marmot::Registry::TransportMaterialFactory::createMaterial( "NOSUCHMATERIAL",
                                                                                           linearIsotropicMaterial,
                                                                                           3,
                                                                                           1 );
                             } ),
                             "an unknown transport material must be rejected" );

    // a solid section on a plane element
    throwExceptionOnFailure( throws( [] {
                               TransportFiniteElement< 2, 4, 1 >( 1,
                                                                  Quadrature::gaussExactToDegree( 2 ),
                                                                  TransportFiniteElement< 2, 4, 1 >::Solid,
                                                                  { "temperature" } );
                             } ),
                             "a solid section on a two-dimensional element must be rejected" );

    // a plane element without thickness
    {
      std::unique_ptr< Element > el( Marmot::Registry::ElementFactory::createElement( "DC2D4T", 1 ) );
      const VectorXd             coords = distortedCoordinates( elementTypes[1], false );
      el->assignNodeCoordinates( coords.data() );
      throwExceptionOnFailure( throws( [&] { el->initializeYourself(); } ),
                               "a plane element without thickness must be rejected" );
    }

    // an axisymmetric element reaching into r < 0
    {
      std::unique_ptr< Element > el( Marmot::Registry::ElementFactory::createElement( "DCAX4T", 1 ) );
      const VectorXd             coords = distortedCoordinates( elementTypes[1], false ); // x in [-1, 1]
      el->assignNodeCoordinates( coords.data() );
      throwExceptionOnFailure( throws( [&] { el->initializeYourself(); } ),
                               "an axisymmetric element at r < 0 must be rejected" );
    }

    // an inverted element (clockwise node order)
    {
      std::unique_ptr< Element > el( Marmot::Registry::ElementFactory::createElement( "DC2D4T", 1 ) );
      const double               coords[] = { 0, 0, 0, 1, 1, 1, 1, 0 };
      std::vector< double >      props{ thickness };
      el->assignNodeCoordinates( coords );
      el->assignProperty( Marmot::Discretizations::ElementProperties( props.data(), 1 ) );
      throwExceptionOnFailure( throws( [&] { el->initializeYourself(); } ), "an inverted element must be rejected" );
    }

    // a field dependent boundary load with the wrong number of nodal values
    {
      const VectorXd coords = distortedCoordinates( elementTypes[1], false );
      const auto     shape  = Marmot::Discretizations::Elements::getElementShapeByMetric( 2, 4 );
      Marmot::Discretizations::Elements::BoundaryElement be( shape, 1, 2, coords, false );
      VectorXd                                           P;
      MatrixXd                                           dP;
      const auto density = []( double, double& value, double& derivative ) { value = derivative = 0.0; };
      throwExceptionOnFailure( throws( [&] {
                                 be.computeFieldDependentScalarLoadVector( VectorXd::Zero( 5 ), density, P, dP );
                               } ),
                               "a boundary field of the wrong size must be rejected" );
    }

    // the defaults of the material base class: storage, and zero initial state variables
    const MinimalTransportMaterial material( linearIsotropicMaterial, 3, 1 );
    throwExceptionOnFailure( material.hasStorage(), "a transport material has storage by default" );
    MinimalTransportMaterial initialized( linearIsotropicMaterial, 3, 1 );
    double                   stateVars[] = { 1.0, 2.0, 3.0 };
    initialized.initializeYourself( stateVars, 3 );
    throwExceptionOnFailure( stateVars[0] == 0.0 && stateVars[1] == 0.0 && stateVars[2] == 0.0,
                             "the default initialization zeroes the state variables" );
    std::cout << "invalid input is rejected, the material defaults hold" << std::endl;
  }

} // namespace

int main()
{
  for ( const auto& type : elementTypes ) {
    const auto element = std::unique_ptr< Element >( Marmot::Registry::ElementFactory::createElement( type.name, 1 ) );
    const auto fields  = element->getNodeFields();
    throwExceptionOnFailure( static_cast< int >( fields.size() ) == type.nNodes && fields[0].size() == 1 &&
                               fields[0][0] == "temperature",
                             type.name + ": one field 'temperature' per node expected" );

    testTangent( type );
    testLoadTangents( type );
    testResultants( type );
    testLinearFieldConsistency( type );
    testStorage( type );
    testInterface( type );
  }
  testInvalidInput();
  return 0;
}
