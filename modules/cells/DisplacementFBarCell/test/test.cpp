/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include "Marmot/DisplacementFBarCell.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  // COMPRESSIBLENEOHOOKE: K, G, rho
  const std::vector< double > matProps = { 3500., 1500., 2.0 };

  // a cell geometry: Lagrangian cells on [0,2] x [0,1] (x [0,1]), optionally distorted, and B-spline cells of order p
  // on the unit knot span [p, p+1]^nDim (control points at the Greville abscissae, see the DisplacementCell test)
  template < int nDim >
  struct CellGeometry {
    std::string                      kind; // "Quad4", "Hexa8", "BSpline/<p>" or "BSpline/3D/<p>"
    std::vector< double >            nodes;
    std::vector< double >            knots;
    Eigen::Matrix< double, nDim, 1 > lower, upper;
    bool                             distorted = false;
  };

  template < int nDim >
  CellGeometry< nDim > lagrangian( bool distorted = false )
  {
    CellGeometry< nDim > g;
    g.distorted = distorted;
    if constexpr ( nDim == 2 ) {
      g.kind  = "Quad4";
      g.nodes = { 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, 1.0 };
      if ( distorted ) // a trapezoid with a skewed top
        g.nodes = { 0.0, 0.0, 2.0, 0.0, 1.7, 1.2, 0.2, 0.9 };
    }
    else {
      g.kind = "Hexa8";
      for ( double z : { 0.0, 1.0 } )
        for ( auto [x, y] : std::vector< std::pair< double, double > >{ { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } } )
          g.nodes.insert( g.nodes.end(), { x, y, z } );
      if ( distorted ) {
        g.nodes[3 * 6 + 0] = 2.3;  // move node 6
        g.nodes[3 * 6 + 2] = 1.2;
        g.nodes[3 * 4 + 1] = -0.1; // and node 4
      }
    }
    Eigen::Map< const Eigen::Matrix< double, nDim, Eigen::Dynamic > > X( g.nodes.data(),
                                                                         nDim,
                                                                         static_cast< int >( g.nodes.size() ) / nDim );
    g.lower = X.rowwise().minCoeff();
    g.upper = X.rowwise().maxCoeff();
    return g;
  }

  template < int nDim >
  CellGeometry< nDim > bSpline( int p )
  {
    CellGeometry< nDim > g;
    g.kind = std::string( "BSpline/" ) + ( nDim == 3 ? "3D/" : "" ) + std::to_string( p );

    const int             nKnots = 2 * p + 2;
    std::vector< double > z( nKnots );
    for ( int k = 0; k < nKnots; k++ )
      z[k] = k;
    for ( int d = 0; d < nDim; d++ )
      g.knots.insert( g.knots.end(), z.begin(), z.end() );

    std::vector< double > greville( p + 1 );
    for ( int i = 0; i <= p; i++ ) {
      greville[i] = 0;
      for ( int k = 1; k <= p; k++ )
        greville[i] += z[i + k] / p;
    }
    const int nN = p + 1;
    for ( int r = 0; r < ( nDim == 3 ? nN : 1 ); r++ )
      for ( int q = 0; q < nN; q++ )
        for ( int pp = 0; pp < nN; pp++ ) {
          g.nodes.push_back( greville[pp] );
          g.nodes.push_back( greville[q] );
          if constexpr ( nDim == 3 )
            g.nodes.push_back( greville[r] );
        }
    g.lower.setConstant( p );
    g.upper.setConstant( p + 1 );
    return g;
  }

  // where the material points sit, in fractions of the bounding box
  enum class Points {
    Symmetric,     ///< 2^nDim points at 1/4 and 3/4
    SingleOffset,  ///< a single point far off the cell center
    ThreeOffCenter ///< three points clustered off the center
  };

  // one cell with its material points, as the host (EdelweissMeshfree) drives it
  template < int nDim >
  struct Setup {
    std::unique_ptr< Marmot::Discretizations::Cells::Cell >                                  cell;
    std::vector< std::unique_ptr< Marmot::Discretizations::MaterialPoints::MaterialPoint > > mps;
    std::vector< std::vector< double > >                                                     stateVars;
    Marmot::Discretizations::MaterialSection section{ "COMPRESSIBLENEOHOOKE", matProps.data(), (int)matProps.size() };
    int                                      nDof;

    Setup( const CellGeometry< nDim >& g, Points points, bool fBar = true )
    {
      const std::string name = std::string( fBar ? "DisplacementFBar/" : "Displacement/" ) + g.kind;
      if ( g.knots.empty() )
        cell.reset(
          Marmot::Registry::CellFactory::createCell( name, 1, g.nodes.data(), static_cast< int >( g.nodes.size() ) ) );
      else
        cell.reset( Marmot::Registry::CellFactory::createBSplineCell( name,
                                                                      1,
                                                                      g.nodes.data(),
                                                                      static_cast< int >( g.nodes.size() ),
                                                                      g.knots.data(),
                                                                      static_cast< int >( g.knots.size() ) ) );
      nDof = cell->getNDofPerCell();

      std::vector< Eigen::Matrix< double, nDim, 1 > > fractions;
      switch ( points ) {
      case Points::Symmetric:
        for ( int k = 0; k < ( 1 << nDim ); k++ ) {
          Eigen::Matrix< double, nDim, 1 > x;
          for ( int d = 0; d < nDim; d++ )
            x[d] = ( k >> d ) & 1 ? 0.75 : 0.25;
          fractions.push_back( x );
        }
        break;
      case Points::SingleOffset:
        fractions.push_back( Eigen::Matrix< double, nDim, 1 >::LinSpaced( 0.85, 0.15 ) );
        break;
      case Points::ThreeOffCenter:
        fractions.push_back( Eigen::Matrix< double, nDim, 1 >::LinSpaced( 0.80, 0.60 ) );
        fractions.push_back( Eigen::Matrix< double, nDim, 1 >::LinSpaced( 0.65, 0.85 ) );
        fractions.push_back( Eigen::Matrix< double, nDim, 1 >::Constant( 0.9 ) );
        break;
      }

      // the points in the (possibly distorted) cell, mapped from the bounding box fractions: for distorted cells the
      // inverse mapping is used through the interpolation vector of the cell, so that the points are inside
      const int nNodes  = cell->getNNodes();
      int       mpLabel = 1;
      for ( const auto& frac : fractions ) {
        Eigen::Matrix< double, nDim, 1 > x = g.lower + frac.cwiseProduct( g.upper - g.lower );
        if ( g.distorted ) { // bilinear/trilinear map of the fractions: x = sum_A N_A(xi) X_A, xi = 2 frac - 1
          Eigen::Matrix< double, nDim, 1 > xi = 2 * frac - Eigen::Matrix< double, nDim, 1 >::Ones();
          x.setZero();
          for ( int A = 0; A < nNodes; A++ ) {
            // the corner signs of the standard node numbering
            const double s[8][3] = { { -1, -1, -1 },
                                     { 1, -1, -1 },
                                     { 1, 1, -1 },
                                     { -1, 1, -1 },
                                     { -1, -1, 1 },
                                     { 1, -1, 1 },
                                     { 1, 1, 1 },
                                     { -1, 1, 1 } };
            double       N       = 1;
            for ( int d = 0; d < nDim; d++ )
              N *= ( 1 + s[A][d] * xi[d] ) / 2;
            x += N * Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( &g.nodes[A * nDim] );
          }
        }
        mps.emplace_back( Marmot::Registry::MaterialPointFactory::createMaterialPoint( nDim == 2
                                                                                         ? "Displacement/PlaneStrain"
                                                                                         : "Displacement/3D",
                                                                                       mpLabel++,
                                                                                       x.data(),
                                                                                       nDim,
                                                                                       0.1 ) );
      }

      for ( auto& mp : mps ) {
        mp->assignMaterial( section );
        stateVars.emplace_back( mp->getNumberOfRequiredStateVars(), 0.0 );
        mp->assignStateVars( stateVars.back().data(), static_cast< int >( stateVars.back().size() ) );
        mp->initializeYourself();
      }

      std::vector< Marmot::Discretizations::MaterialPoints::MaterialPoint* > raw;
      for ( auto& mp : mps )
        raw.push_back( mp.get() );
      cell->assignMaterialPoints( raw );
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& dQ )
    {
      for ( auto& mp : mps )
        mp->prepareYourself( 1.0, 1.0 );
      cell->interpolateFieldsToMaterialPoints( dQ.data() );
      for ( auto& mp : mps )
        mp->computeYourself( 1.0, 1.0 );

      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      cell->computeMaterialPointKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
      return { P, K };
    }

    // evaluate without committing: the material history is restored afterwards
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > trial( const Eigen::VectorXd& dQ )
    {
      const auto backup = stateVars;
      auto       result = kernels( dQ );
      stateVars         = backup;
      for ( size_t i = 0; i < mps.size(); i++ )
        mps[i]->assignStateVars( stateVars[i].data(), static_cast< int >( stateVars[i].size() ) );
      return result;
    }

    void accept( const Eigen::VectorXd& dQ )
    {
      kernels( dQ );
      for ( auto& mp : mps )
        mp->acceptStateAndPosition();
    }
  };

  Eigen::VectorXd increment( int nDof, double scale )
  {
    Eigen::VectorXd dQ( nDof );
    for ( int i = 0; i < nDof; i++ )
      dQ[i] = scale * std::sin( 1.7 * ( i + 1 ) );
    return dQ;
  }

  // a homogeneous increment dq_A = H X_A
  template < int nDim >
  Eigen::VectorXd homogeneousIncrement( const CellGeometry< nDim >& g )
  {
    Eigen::Matrix< double, nDim, nDim > H;
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        H( i, j ) = 0.03 * std::cos( 1.3 * ( 1 + i + nDim * j ) );
    const int       nNodes = static_cast< int >( g.nodes.size() ) / nDim;
    Eigen::VectorXd dQ( nDim * nNodes );
    for ( int A = 0; A < nNodes; A++ )
      dQ.segment< nDim >( nDim * A ) = H * Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( &g.nodes[A * nDim] );
    return dQ;
  }

  // relative error of K against central differences of f along directions D (the identity up to 64 dofs, three
  // random unit directions beyond, see the DisplacementCell test)
  double tangentError( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                       const Eigen::VectorXd&                                            Q0,
                       const Eigen::MatrixXd&                                            K )
  {
    const int       n = static_cast< int >( Q0.size() );
    Eigen::MatrixXd D;
    if ( n <= 64 )
      D = Eigen::MatrixXd::Identity( n, n );
    else {
      std::mt19937                             rng( 1 );
      std::uniform_real_distribution< double > uniform( -1.0, 1.0 );
      D.resize( n, 3 );
      for ( int k = 0; k < 3; k++ ) {
        for ( int j = 0; j < n; j++ )
          D( j, k ) = uniform( rng );
        D.col( k ).normalize();
      }
    }
    const double    h = 1e-7;
    Eigen::MatrixXd KD( n, D.cols() );
    for ( int j = 0; j < D.cols(); j++ )
      KD.col( j ) = ( f( Q0 + h * D.col( j ) ) - f( Q0 - h * D.col( j ) ) ) / ( 2 * h );
    return ( K * D - KD ).norm() / KD.norm();
  }

  std::string label( const std::string& kind, Points points, bool distorted )
  {
    const char* p[] = { "symmetric points", "single off-center point", "three off-center points" };
    return Marmot::Utilities::MakeString() << "DisplacementFBar/" << kind << ( distorted ? " (distorted)" : "" ) << ", "
                                           << p[static_cast< int >( points )];
  }

  // the exact tangent, from the undeformed state and after an accepted finite increment
  template < int nDim >
  void checkTangent( const CellGeometry< nDim >& g, Points points )
  {
    for ( bool afterStep : { false, true } ) {
      Setup< nDim > s( g, points );
      if ( afterStep )
        s.accept( increment( s.nDof, 0.04 ) );
      const Eigen::VectorXd dQ = increment( s.nDof, afterStep ? 0.02 : 0.05 );
      const auto [P, K]        = s.trial( dQ );
      const double err         = tangentError( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ, K );
      throwExceptionOnFailure( err < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << label( g.kind, points, g.distorted ) << ( afterStep ? ", second" : ", first" )
                                 << " step: tangent inconsistent, relative error " << err );
    }
  }

  // the F-bar projection: all points share the volume change of the center, and the cell differs from the
  // DisplacementCell for a nonhomogeneous increment
  template < int nDim >
  void checkProjection( const CellGeometry< nDim >& g )
  {
    Setup< nDim >         fbar( g, Points::ThreeOffCenter );
    Setup< nDim >         plain( g, Points::ThreeOffCenter, false );
    const Eigen::VectorXd dQ = increment( fbar.nDof, 0.05 );
    const auto            Pf = fbar.trial( dQ ).first;
    const auto            Pp = plain.trial( dQ ).first;
    fbar.kernels( dQ );

    std::vector< double > dets;
    for ( auto& mp : fbar.mps ) {
      Eigen::Map< const Eigen::Matrix3d > dF( mp->getStateView( "delta deformation gradient" ).stateLocation );
      dets.push_back( dF.determinant() );
    }
    for ( double d : dets )
      throwExceptionOnFailure( std::abs( d - dets[0] ) < 1e-13 && std::abs( d - 1 ) > 1e-3,
                               label( g.kind, Points::ThreeOffCenter, g.distorted ) +
                                 ": the points must share the volume change of the center" );
    throwExceptionOnFailure( ( Pf - Pp ).norm() > 1e-6 * Pp.norm(),
                             label( g.kind, Points::ThreeOffCenter, g.distorted ) +
                               ": F-bar must change the response of a nonhomogeneous increment" );
  }

  // reduction to the DisplacementCell: a single point (residual and tangent, also after a step), and a homogeneous
  // increment (residual)
  template < int nDim >
  void checkReduction( const CellGeometry< nDim >& g )
  {
    {
      Setup< nDim > fbar( g, Points::SingleOffset );
      Setup< nDim > plain( g, Points::SingleOffset, false );
      for ( bool afterStep : { false, true } ) {
        if ( afterStep ) {
          fbar.accept( increment( fbar.nDof, 0.04 ) );
          plain.accept( increment( plain.nDof, 0.04 ) );
        }
        const Eigen::VectorXd dQ = increment( fbar.nDof, 0.03 );
        const auto [Pf, Kf]      = fbar.trial( dQ );
        const auto [Pp, Kp]      = plain.trial( dQ );
        throwExceptionOnFailure( ( Pf - Pp ).norm() <= 1e-12 * Pp.norm() && ( Kf - Kp ).norm() <= 1e-12 * Kp.norm(),
                                 label( g.kind, Points::SingleOffset, g.distorted ) +
                                   ": a single point must reproduce the DisplacementCell" );
      }
    }
    for ( Points points : { Points::Symmetric, Points::ThreeOffCenter } ) {
      Setup< nDim >         fbar( g, points );
      Setup< nDim >         plain( g, points, false );
      const Eigen::VectorXd dQ = homogeneousIncrement( g );
      const auto            Pf = fbar.trial( dQ ).first;
      const auto            Pp = plain.trial( dQ ).first;
      throwExceptionOnFailure( Pp.norm() > 0 && ( Pf - Pp ).norm() <= 1e-10 * Pp.norm(),
                               label( g.kind, points, g.distorted ) +
                                 ": a homogeneous increment must reproduce the DisplacementCell" );
    }
  }

  // a finite rigid rotation does not load the cell; the follower pressure has an exact tangent
  template < int nDim >
  void checkRotationAndPressure( const CellGeometry< nDim >& g )
  {
    const int nNodes = static_cast< int >( g.nodes.size() ) / nDim;
    {
      Setup< nDim >   s( g, Points::ThreeOffCenter );
      const double    phi = 0.5;
      Eigen::Matrix3d R   = Eigen::Matrix3d::Identity();
      R.topLeftCorner< 2, 2 >() << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );
      Eigen::VectorXd dQ( s.nDof );
      for ( int A = 0; A < nNodes; A++ ) {
        Eigen::Vector3d X              = Eigen::Vector3d::Zero();
        X.head< nDim >()               = Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( &g.nodes[A * nDim] );
        dQ.segment< nDim >( nDim * A ) = ( ( R - Eigen::Matrix3d::Identity() ) * X ).template head< nDim >();
      }
      const auto P = s.trial( dQ ).first;
      throwExceptionOnFailure( P.norm() < 1e-8,
                               Marmot::Utilities::MakeString()
                                 << label( g.kind, Points::ThreeOffCenter, g.distorted )
                                 << ": rigid rotation must not load the cell, |P| = " << P.norm() );
    }
    {
      Setup< nDim >                          s( g, Points::ThreeOffCenter );
      const Eigen::VectorXd                  dQ       = increment( s.nDof, 0.05 );
      const Eigen::Matrix< double, nDim, 1 > p        = -2.0 * Eigen::Matrix< double, nDim, 1 >::Unit( nDim - 1 );
      const int                              pressure = s.cell->getSupportedDistributedLoadTypes().at( "PRESSURE" );

      auto load = [&]( const Eigen::VectorXd& q ) {
        s.trial( q );
        for ( auto& mp : s.mps )
          mp->prepareYourself( 1.0, 1.0 );
        s.cell->interpolateFieldsToMaterialPoints( q.data() );
        Eigen::VectorXd fExt = Eigen::VectorXd::Zero( s.nDof );
        Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( s.nDof, s.nDof );
        s.cell->computeDistributedLoad( pressure, 1, 2, p.data(), fExt.data(), K.data(), 1.0, 1.0 );
        return std::make_pair( fExt, K );
      };

      const auto [f0, K0] = load( dQ );
      const double err    = tangentError( [&]( const Eigen::VectorXd& q ) { return load( q ).first; }, dQ, K0 );
      throwExceptionOnFailure( f0.norm() > 0 && err < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << label( g.kind, Points::ThreeOffCenter, g.distorted )
                                 << ": pressure load tangent inconsistent, error " << err );
    }
  }

  template < int nDim >
  void checkCell( const CellGeometry< nDim >& g )
  {
    for ( Points points : { Points::Symmetric, Points::SingleOffset, Points::ThreeOffCenter } )
      checkTangent( g, points );
    checkProjection( g );
    checkReduction( g );
    checkRotationAndPressure( g );
  }

  // the projection point (aggregated MPM): at the center of the points it changes nothing; two cells carrying one
  // bilinear field (the second one its extension) and evaluating the volume change at one common point give all their
  // points the same volume change; the tangent stays exact for a point outside the cell; B-spline cells decline
  void checkProjectionPoint()
  {
    const auto g = lagrangian< 2 >();
    {
      Setup< 2 >      a( g, Points::ThreeOffCenter );
      Setup< 2 >      b( g, Points::ThreeOffCenter );
      Eigen::Vector2d center = Eigen::Vector2d::Zero();
      for ( auto& mp : b.mps ) {
        Eigen::Vector2d x;
        mp->getCoordinatesAtCenter( x.data() );
        center += x / b.mps.size();
      }
      throwExceptionOnFailure( b.cell->setProjectionPoint( center.data() ),
                               "DisplacementFBar/Quad4 must use a projection point" );
      const Eigen::VectorXd dQ = increment( a.nDof, 0.05 );
      const auto [Pa, Ka]      = a.trial( dQ );
      const auto [Pb, Kb]      = b.trial( dQ );
      throwExceptionOnFailure( ( Pa - Pb ).norm() <= 1e-12 * Pa.norm() && ( Ka - Kb ).norm() <= 1e-12 * Ka.norm(),
                               "DisplacementFBar/Quad4: the projection point at the center must change nothing" );
    }
    {
      // cell 1 on [0,2] x [0,1], cell 2 on [2,4] x [0,1]; the field u = H x + c xy, common point in cell 1
      CellGeometry< 2 > g2 = g;
      g2.nodes             = { 2.0, 0.0, 4.0, 0.0, 4.0, 1.0, 2.0, 1.0 };
      g2.lower << 2.0, 0.0;
      g2.upper << 4.0, 1.0;
      Setup< 2 >      s1( g, Points::ThreeOffCenter );
      Setup< 2 >      s2( g2, Points::Symmetric );
      Eigen::Matrix2d H;
      H << 0.02, -0.03, 0.01, 0.04;
      const Eigen::Vector2d c( 0.03, -0.02 );
      auto                  field = [&]( const CellGeometry< 2 >& gg ) {
        Eigen::VectorXd dQ( 8 );
        for ( int A = 0; A < 4; A++ ) {
          const Eigen::Vector2d X( gg.nodes[2 * A], gg.nodes[2 * A + 1] );
          dQ.segment< 2 >( 2 * A ) = H * X + c * X[0] * X[1];
        }
        return dQ;
      };
      const Eigen::Vector2d point( 0.7, 0.4 );
      s1.cell->setProjectionPoint( point.data() );
      s2.cell->setProjectionPoint( point.data() );
      s1.kernels( field( g ) );
      s2.kernels( field( g2 ) );
      std::vector< double > dets;
      for ( auto* s : { &s1, &s2 } )
        for ( auto& mp : s->mps ) {
          Eigen::Map< const Eigen::Matrix3d > dF( mp->getStateView( "delta deformation gradient" ).stateLocation );
          dets.push_back( dF.determinant() );
        }
      for ( double d : dets )
        throwExceptionOnFailure( std::abs( d - dets[0] ) < 1e-13 && std::abs( d - 1 ) > 1e-3,
                                 "DisplacementFBar/Quad4: two cells of one field must share the volume change of the "
                                 "common projection point" );

      // exact tangent with the point outside the cell, also after a step
      for ( bool afterStep : { false, true } ) {
        Setup< 2 > s( g2, Points::ThreeOffCenter );
        s.cell->setProjectionPoint( point.data() );
        if ( afterStep ) {
          s.accept( increment( s.nDof, 0.04 ) );
          s.cell->setProjectionPoint( point.data() );
        }
        const Eigen::VectorXd dQ = increment( s.nDof, 0.03 );
        const auto [P, K]        = s.trial( dQ );
        const double err = tangentError( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ, K );
        throwExceptionOnFailure( err < 1e-7,
                                 Marmot::Utilities::MakeString()
                                   << "DisplacementFBar/Quad4, projection point outside the cell: tangent "
                                      "inconsistent, relative error "
                                   << err );
      }
    }
    {
      Setup< 2 >   s( bSpline< 2 >( 2 ), Points::Symmetric );
      const double x[2] = { 0.0, 0.0 };
      throwExceptionOnFailure( !s.cell->setProjectionPoint( x ),
                               "B-spline F-bar cells must decline the projection point" );
    }
  }

  // ---------------------------------------------------------------------------------------------------------------
  // F-bar patch: three Quad4 cells sharing nodes (or one), driven in the host's two passes; global dofs

  // FINITESTRAINJ2PLASTICITY: K, G, fy, fyInf, eta, H, implementationType, density
  const std::vector< double > j2Props = { 833., 385., 6., 6., 0., 0., 1., 1.0 };

  struct Patch {
    // global nodes (0,0) (2,0) (4,0) (0,1) (2,1) (4,1) (0,2) (2,2); cells [0,1,4,3] [1,2,5,4] [3,4,7,6]
    std::vector< double >                                                                    X;
    std::vector< std::array< int, 4 > >                                                      conn;
    std::array< std::vector< double >, 3 >                                                   cellNodes;
    std::vector< std::unique_ptr< Marmot::Discretizations::Cells::Cell > >                   cells;
    std::vector< std::vector< Marmot::Discretizations::MaterialPoints::MaterialPoint* > >    cellMps;
    std::vector< std::unique_ptr< Marmot::Discretizations::MaterialPoints::MaterialPoint > > mps;
    std::vector< std::vector< double > >                                                     stateVars;
    std::unique_ptr< Marmot::Discretizations::MaterialSection >                              section;
    int                                                                                      nDof;
    bool                                                                                     patch;

    Patch( int nCells, bool plastic, bool patchMode = true, bool fBar = true ) : patch( patchMode )
    {
      X    = { 0, 0, 2, 0, 4, 0, 0, 1, 2, 1, 4, 1, 0, 2, 2, 2 };
      conn = { { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 3, 4, 7, 6 } };
      conn.resize( nCells );
      nDof = 16;
      section.reset( plastic ? new Marmot::Discretizations::MaterialSection( "FINITESTRAINJ2PLASTICITY",
                                                                             j2Props.data(),
                                                                             (int)j2Props.size() )
                             : new Marmot::Discretizations::MaterialSection( "COMPRESSIBLENEOHOOKE",
                                                                             matProps.data(),
                                                                             (int)matProps.size() ) );
      // material points in fractions of each cell, off the cell centers
      const std::vector< std::vector< std::array< double, 2 > > > fractions =
        { { { 0.80, 0.60 }, { 0.65, 0.85 }, { 0.9, 0.9 } },
          { { 0.2, 0.3 }, { 0.7, 0.2 }, { 0.3, 0.75 }, { 0.8, 0.8 } },
          { { 0.15, 0.85 }, { 0.3, 0.6 } } };
      int label_ = 1;
      for ( int c = 0; c < nCells; c++ ) {
        auto& nodes = cellNodes[c]; // the cell keeps a view of its node coordinates
        for ( int a : conn[c] )
          nodes.insert( nodes.end(), { X[2 * a], X[2 * a + 1] } );
        cells.emplace_back(
          Marmot::Registry::CellFactory::createCell( fBar ? "DisplacementFBar/Quad4" : "Displacement/Quad4",
                                                     c + 1,
                                                     nodes.data(),
                                                     8 ) );
        cellMps.emplace_back();
        for ( const auto& f : fractions[c] ) {
          const double x[2] = { nodes[0] + f[0] * 2.0, nodes[1] + f[1] * 1.0 };
          mps.emplace_back( Marmot::Registry::MaterialPointFactory::createMaterialPoint( "Displacement/PlaneStrain",
                                                                                         label_++,
                                                                                         x,
                                                                                         2,
                                                                                         0.1 + 0.05 * c ) );
          cellMps.back().push_back( mps.back().get() );
        }
      }
      for ( auto& mp : mps ) {
        mp->assignMaterial( *section );
        stateVars.emplace_back( mp->getNumberOfRequiredStateVars(), 0.0 );
        mp->assignStateVars( stateVars.back().data(), static_cast< int >( stateVars.back().size() ) );
        mp->initializeYourself();
      }
      assign();
    }

    void assign()
    {
      for ( size_t c = 0; c < cells.size(); c++ )
        cells[c]->assignMaterialPoints( cellMps[c] );
    }

    Eigen::VectorXd gather( const Eigen::VectorXd& U, size_t c ) const
    {
      Eigen::VectorXd q( 8 );
      for ( int a = 0; a < 4; a++ )
        q.segment< 2 >( 2 * a ) = U.segment< 2 >( 2 * conn[c][a] );
      return q;
    }

    // the host's evaluation: pass 1 (patch value), interpolation, material points, kernels, patch coupling
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& U, double* JA = nullptr )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof ), dv = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
      if ( patch ) {
        double v = 0, w = 0;
        for ( size_t c = 0; c < cells.size(); c++ ) {
          double          vc, wc;
          Eigen::VectorXd dvc( 8 );
          cells[c]->computeProjectionContribution( gather( U, c ).data(), &vc, &wc, dvc.data() );
          v += vc;
          w += wc;
          for ( int a = 0; a < 4; a++ )
            dv.segment< 2 >( 2 * conn[c][a] ) += dvc.segment< 2 >( 2 * a );
        }
        for ( auto& cell : cells )
          cell->setProjectedValue( v / w );
        dv /= w;
        if ( JA )
          *JA = v / w;
      }
      for ( auto& mp : mps )
        mp->prepareYourself( 1.0, 1.0 );
      for ( size_t c = 0; c < cells.size(); c++ )
        cells[c]->interpolateFieldsToMaterialPoints( gather( U, c ).data() );
      for ( auto& mp : mps )
        mp->computeYourself( 1.0, 1.0 );
      for ( size_t c = 0; c < cells.size(); c++ ) {
        Eigen::VectorXd Pc = Eigen::VectorXd::Zero( 8 ), a( 8 );
        Eigen::MatrixXd Kc = Eigen::MatrixXd::Zero( 8, 8 );
        const auto      q  = gather( U, c );
        cells[c]->computeMaterialPointKernels( q.data(), Pc.data(), Kc.data(), 1.0, 1.0 );
        Eigen::VectorXd aG = Eigen::VectorXd::Zero( nDof );
        if ( patch )
          cells[c]->getProjectionSensitivity( a.data() );
        for ( int A = 0; A < 4; A++ ) {
          P.segment< 2 >( 2 * conn[c][A] ) += Pc.segment< 2 >( 2 * A );
          if ( patch )
            aG.segment< 2 >( 2 * conn[c][A] ) += a.segment< 2 >( 2 * A );
          for ( int B = 0; B < 4; B++ )
            K.block< 2, 2 >( 2 * conn[c][A], 2 * conn[c][B] ) += Kc.block< 2, 2 >( 2 * A, 2 * B );
        }
        if ( patch )
          K += aG * dv.transpose();
      }
      return { P, K };
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > trial( const Eigen::VectorXd& U, double* JA = nullptr )
    {
      const auto backup = stateVars;
      auto       result = kernels( U, JA );
      stateVars         = backup;
      for ( size_t i = 0; i < mps.size(); i++ )
        mps[i]->assignStateVars( stateVars[i].data(), static_cast< int >( stateVars[i].size() ) );
      return result;
    }

    void accept( const Eigen::VectorXd& U )
    {
      kernels( U );
      for ( auto& mp : mps )
        mp->acceptStateAndPosition();
      assign(); // a new increment: the host reassigns (the points stay in their cells)
    }

    double volumeChange( size_t k ) const
    {
      Eigen::Map< const Eigen::Matrix3d > dF( mps[k]->getStateView( "delta deformation gradient" ).stateLocation );
      return dF.determinant();
    }
  };

  void checkPatch()
  {
    for ( int nCells : { 1, 3 } )
      for ( bool plastic : { false, true } )
        for ( bool afterStep : { false, true } ) {
          const std::string what = Marmot::Utilities::MakeString()
                                   << "DisplacementFBar/Quad4 patch of " << nCells << " cell(s), "
                                   << ( plastic ? "J2 plastic" : "elastic" ) << ( afterStep ? ", second" : ", first" )
                                   << " step";
          Patch        s( nCells, plastic );
          const double scale = plastic ? 0.25 : 1.0;
          if ( afterStep )
            s.accept( increment( s.nDof, 0.03 * scale ) );
          const Eigen::VectorXd U = increment( s.nDof, ( afterStep ? 0.02 : 0.04 ) * scale );
          double                JA;
          const auto [P, K]            = s.trial( U, &JA );
          const auto            backup = s.stateVars;
          std::vector< double > alphaPn;
          if ( plastic )
            for ( auto& mp : s.mps )
              alphaPn.push_back( mp->getStateView( "alphaP" ).stateLocation[0] );
          s.kernels( U ); // the state of the trial, to read the volume changes
          double v = 0, w = 0;
          for ( size_t k = 0; k < s.mps.size(); k++ ) {
            throwExceptionOnFailure( std::abs( s.volumeChange( k ) - JA ) < 1e-13 && std::abs( JA - 1 ) > 1e-4,
                                     Marmot::Utilities::MakeString()
                                       << what << ": every point must receive the volume change of the patch, "
                                       << s.volumeChange( k ) << " vs " << JA );
            Eigen::Map< const Eigen::Matrix3d > Fn( s.mps[k]->getStateView( "deformation gradient" ).stateLocation );
            v += s.mps[k]->getVolumeUndeformed() * Fn.determinant();
            w += s.mps[k]->getVolumeUndeformed() * Fn.determinant() / JA;
          }
          // the patch keeps its volume: sum V_n dJ_A = sum V_n det dF_p (checked through the mean of the J's)
          throwExceptionOnFailure( std::abs( v / w - JA ) < 1e-12, what + ": patch volume" );
          if ( plastic ) {
            bool yielded = false;
            for ( size_t k = 0; k < s.mps.size(); k++ )
              yielded |= s.mps[k]->getStateView( "alphaP" ).stateLocation[0] > alphaPn[k] + 1e-8;
            throwExceptionOnFailure( yielded, what + ": the test must reach the plastic regime" );
          }
          s.stateVars = backup;
          for ( size_t i = 0; i < s.mps.size(); i++ )
            s.mps[i]->assignStateVars( s.stateVars[i].data(), static_cast< int >( s.stateVars[i].size() ) );
          const double err = tangentError( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, U, K );
          throwExceptionOnFailure( err < 1e-7,
                                   Marmot::Utilities::MakeString()
                                     << what << ": tangent of the patch inconsistent, relative error " << err );
        }

    // at a zero increment from the reference state the tangent of the patch is symmetric positive semidefinite
    {
      Patch                 s( 3, false );
      const Eigen::VectorXd U                               = Eigen::VectorXd::Zero( s.nDof );
      const auto [P, K]                                     = s.trial( U );
      const double                                     asym = ( K - K.transpose() ).norm() / K.norm();
      Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eig( 0.5 * ( K + K.transpose() ) );
      throwExceptionOnFailure( asym < 1e-12 && eig.eigenvalues().minCoeff() > -1e-9 * eig.eigenvalues().maxCoeff(),
                               Marmot::Utilities::MakeString()
                                 << "F-bar patch at zero increment: asymmetry " << asym << ", lowest eigenvalue "
                                 << eig.eigenvalues().minCoeff() );
    }

    // rigid rotation: no load; homogeneous increment: the DisplacementCell
    {
      Patch           s( 3, false );
      const double    phi = 0.5;
      Eigen::Matrix2d R;
      R << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );
      Eigen::Matrix2d H;
      H << 0.02, -0.03, 0.01, 0.04;
      Eigen::VectorXd Ur( s.nDof ), Uh( s.nDof );
      for ( int a = 0; a < 8; a++ ) {
        const Eigen::Vector2d Xa( s.X[2 * a], s.X[2 * a + 1] );
        Ur.segment< 2 >( 2 * a ) = ( R - Eigen::Matrix2d::Identity() ) * Xa;
        Uh.segment< 2 >( 2 * a ) = H * Xa;
      }
      const auto Pr = s.trial( Ur ).first;
      throwExceptionOnFailure( Pr.norm() < 1e-8,
                               Marmot::Utilities::MakeString()
                                 << "F-bar patch: rigid rotation must not load the patch, |P| = " << Pr.norm() );
      Patch      plain( 3, false, false, false );
      const auto Ph = s.trial( Uh ).first;
      const auto Pp = plain.trial( Uh ).first;
      throwExceptionOnFailure( ( Ph - Pp ).norm() <= 1e-10 * Pp.norm(),
                               "F-bar patch: a homogeneous increment must reproduce the DisplacementCell" );
    }

    // a patch of one cell is the cell-average variant: the volume-weighted mean of the det dF_p of the cell
    {
      Patch                 s( 1, false );
      Patch                 plain( 1, false, false, false );
      const Eigen::VectorXd U = increment( s.nDof, 0.04 );
      double                JA;
      s.trial( U, &JA );
      plain.kernels( U );
      double v = 0, w = 0;
      for ( size_t k = 0; k < plain.mps.size(); k++ ) {
        v += plain.mps[k]->getVolumeUndeformed() * plain.volumeChange( k );
        w += plain.mps[k]->getVolumeUndeformed();
      }
      throwExceptionOnFailure( std::abs( JA - v / w ) < 1e-13,
                               "F-bar patch of one cell: the volume change must be the mean over the cell" );
    }

    // the B-spline cells support the patch too
    {
      Setup< 2 > s( bSpline< 2 >( 2 ), Points::ThreeOffCenter );
      throwExceptionOnFailure( s.cell->hasPatchProjection(), "B-spline F-bar cells must support the patch" );
    }
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > testFunctions = { []() { checkCell< 2 >( lagrangian< 2 >() ); },
                                                           []() { checkCell< 2 >( lagrangian< 2 >( true ) ); },
                                                           []() { checkCell< 3 >( lagrangian< 3 >() ); },
                                                           []() { checkCell< 3 >( lagrangian< 3 >( true ) ); },
                                                           []() { checkProjectionPoint(); },
                                                           []() { checkPatch(); } };
  for ( int p = 1; p <= 3; p++ ) {
    testFunctions.push_back( [p]() { checkCell< 2 >( bSpline< 2 >( p ) ); } );
    testFunctions.push_back( [p]() { checkCell< 3 >( bSpline< 3 >( p ) ); } );
  }
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
