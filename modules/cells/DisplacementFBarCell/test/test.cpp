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

} // namespace

int main()
{
  std::vector< std::function< void() > > testFunctions = { []() { checkCell< 2 >( lagrangian< 2 >() ); },
                                                           []() { checkCell< 2 >( lagrangian< 2 >( true ) ); },
                                                           []() { checkCell< 3 >( lagrangian< 3 >() ); },
                                                           []() { checkCell< 3 >( lagrangian< 3 >( true ) ); },
                                                           []() { checkProjectionPoint(); } };
  for ( int p = 1; p <= 3; p++ ) {
    testFunctions.push_back( [p]() { checkCell< 2 >( bSpline< 2 >( p ) ); } );
    testFunctions.push_back( [p]() { checkCell< 3 >( bSpline< 3 >( p ) ); } );
  }
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
