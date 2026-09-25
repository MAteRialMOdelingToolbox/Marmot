#include "Marmot/GradientEnhancedFiniteStrainCell.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };
  constexpr int               iRho     = 5;

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

  // a cell and its geometry: Lagrangian box cells on [0,2] x [0,1] (x [0,1]), and B-spline cells of order p on the
  // unit knot span [p, p+1]^nDim of the uniform knot vector 0, 1, ..., 2p+1, with the control points at the
  // Greville abscissae (which gives the B-spline cell linear precision)
  template < int nDim >
  struct CellGeometry {
    std::string                      name;
    std::vector< double >            nodes;
    std::vector< double >            knots; // empty for Lagrangian cells
    Eigen::Matrix< double, nDim, 1 > lower, upper;
  };

  template < int nDim >
  CellGeometry< nDim > lagrangian()
  {
    CellGeometry< nDim > g;
    if constexpr ( nDim == 2 ) {
      g.name  = "GradientEnhancedFiniteStrain/Quad4";
      g.nodes = { 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, 1.0 };
    }
    else {
      g.name = "GradientEnhancedFiniteStrain/Hexa8";
      for ( double z : { 0.0, 1.0 } )
        for ( auto [x, y] : std::vector< std::pair< double, double > >{ { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } } )
          g.nodes.insert( g.nodes.end(), { x, y, z } );
    }
    g.lower.setZero();
    g.upper.setOnes();
    g.upper[0] = 2.0;
    return g;
  }

  template < int nDim >
  CellGeometry< nDim > bSpline( int p )
  {
    CellGeometry< nDim > g;
    g.name = std::string( "GradientEnhancedFiniteStrain/BSpline/" ) + ( nDim == 3 ? "3D/" : "" ) + std::to_string( p );

    const int             nKnots = 2 * p + 2;
    std::vector< double > z( nKnots );
    for ( int k = 0; k < nKnots; k++ )
      z[k] = k;
    for ( int d = 0; d < nDim; d++ ) // column major: one knot vector per direction
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
          if ( nDim == 3 )
            g.nodes.push_back( greville[r] );
        }
    g.lower.setConstant( p );
    g.upper.setConstant( p + 1 );
    return g;
  }

  // one cell with 2^nDim material points, as the host (EdelweissMeshfree) drives it
  template < int nDim >
  struct Setup {
    std::unique_ptr< MarmotCell >                         cell;
    std::vector< std::unique_ptr< MarmotMaterialPoint > > mps;
    std::vector< std::vector< double > >                  stateVars;
    MarmotMaterialSection                                 section{ "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                   matProps.data(),
                                   (int)matProps.size() };
    int                                                   nDof;
    double                                                mpVolume;

    Setup( const CellGeometry< nDim >& g )
    {
      if ( g.knots.empty() )
        cell.reset( MarmotLibrary::MarmotCellFactory::createCell( g.name, 1, g.nodes.data(), g.nodes.size() ) );
      else
        cell.reset( MarmotLibrary::MarmotCellFactory::createBSplineCell( g.name,
                                                                         1,
                                                                         g.nodes.data(),
                                                                         g.nodes.size(),
                                                                         g.knots.data(),
                                                                         g.knots.size() ) );
      nDof     = cell->getNDofPerCell();
      mpVolume = ( g.upper - g.lower ).prod() / std::pow( 2, nDim );

      int label = 1;
      for ( int k = 0; k < ( 1 << nDim ); k++ ) {
        Eigen::Matrix< double, nDim, 1 > x;
        for ( int d = 0; d < nDim; d++ )
          x[d] = g.lower[d] + ( ( k >> d ) & 1 ? 0.75 : 0.25 ) * ( g.upper[d] - g.lower[d] );
        mps.emplace_back(
          MarmotLibrary::MarmotMaterialPointFactory::createMaterialPoint( nDim == 2
                                                                            ? "GradientEnhancedFiniteStrain/PlaneStrain"
                                                                            : "GradientEnhancedFiniteStrain/3D",
                                                                          label++,
                                                                          x.data(),
                                                                          nDim,
                                                                          mpVolume ) );
      }

      for ( auto& mp : mps ) {
        mp->assignMaterial( section );
        stateVars.emplace_back( mp->getNumberOfRequiredStateVars(), 0.0 );
        mp->assignStateVars( stateVars.back().data(), stateVars.back().size() );
        mp->initializeYourself();
      }

      std::vector< MarmotMaterialPoint* > raw;
      for ( auto& mp : mps )
        raw.push_back( mp.get() );
      cell->assignMaterialPoints( raw );
    }

    // one MPM increment: reset the increment, interpolate, compute the points, assemble
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
        mps[i]->assignStateVars( stateVars[i].data(), stateVars[i].size() );
      return result;
    }

    void accept( const Eigen::VectorXd& dQ )
    {
      kernels( dQ );
      for ( auto& mp : mps )
        mp->acceptStateAndPosition();
    }
  };

  // blocked cell dof vector [u | N], nU = nDim nNodes
  Eigen::VectorXd increment( int nDof, int nU, double scaleU, double N0 = 0.0 )
  {
    Eigen::VectorXd dQ( nDof );
    for ( int i = 0; i < nU; i++ )
      dQ[i] = scaleU * std::sin( 1.7 * ( i + 1 ) );
    for ( int A = 0; A < nDof - nU; A++ )
      dQ[nU + A] = N0 + 0.3 * N0 * std::cos( 0.9 * ( A + 1 ) );
    return dQ;
  }

  Eigen::MatrixXd numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                    const Eigen::VectorXd&                                            Q0,
                                    int                                                               nU )
  {
    Eigen::MatrixXd numK( f( Q0 ).size(), Q0.size() );
    for ( int j = 0; j < Q0.size(); j++ ) {
      const double    h  = j < nU ? 1e-7 : 1e-9; // the nonlocal field is small
      Eigen::VectorXd Qp = Q0, Qm = Q0;
      Qp[j] += h;
      Qm[j] -= h;
      numK.col( j ) = ( f( Qp ) - f( Qm ) ) / ( 2 * h );
    }
    return numK;
  }

  template < int nDim >
  void checkCell( const CellGeometry< nDim >& g )
  {
    using Vec          = Eigen::Matrix< double, nDim, 1 >;
    const int  nNodes  = g.nodes.size() / nDim;
    const int  nU      = nDim * nNodes;
    const auto nodeAt  = [&]( int A ) { return Eigen::Map< const Vec >( &g.nodes[A * nDim] ); };
    const Vec  inside  = g.lower + 0.37 * ( g.upper - g.lower );
    const Vec  outside = g.upper + 0.1 * Vec::Ones();

    // layout and geometry
    {
      Setup< nDim > s( g );
      throwExceptionOnFailure( s.cell->getNNodes() == nNodes && s.nDof == ( nDim + 1 ) * nNodes, g.name + ": layout" );
      for ( const auto& f : s.cell->getNodeFields() )
        throwExceptionOnFailure( f == std::vector< std::string >{ "displacement", "nonlocal damage" },
                                 g.name + ": node fields" );
      throwExceptionOnFailure( int( s.cell->getDofIndicesPermutationPattern().size() ) == s.nDof,
                               g.name + ": permutation pattern" );
      // (the shape names cover the linear and quadratic cells; the cubic B-spline cells have none)
      const bool cubic = !g.knots.empty() && g.name.back() == '3';
      throwExceptionOnFailure( cubic || !s.cell->getCellShape().empty(), g.name + ": shape" );

      throwExceptionOnFailure( s.cell->isCoordinateInCell( inside.data() ) &&
                                 !s.cell->isCoordinateInCell( outside.data() ),
                               g.name + ": coordinate in cell" );
      Vec lo, hi;
      s.cell->getBoundingBox( lo.data(), hi.data() );
      throwExceptionOnFailure( ( lo - g.lower ).norm() < 1e-14 && ( hi - g.upper ).norm() < 1e-14,
                               g.name + ": bounding box" );

      // partition of unity and linear precision of the interpolation
      Eigen::VectorXd N( nNodes );
      s.cell->getInterpolationVector( N.data(), inside.data() );
      Vec reproduced = Vec::Zero();
      for ( int A = 0; A < nNodes; A++ )
        reproduced += N[A] * nodeAt( A );
      throwExceptionOnFailure( std::abs( N.sum() - 1 ) < 1e-12 && ( reproduced - inside ).norm() < 1e-12,
                               g.name + ": interpolation" );

      const double zero[3] = { 0, 0, 0 };
      double       dummy[1];
      throwExceptionOnFailure( throws( [&]() { s.cell->computeBodyLoad( -1, zero, dummy, dummy, 0., 1. ); } ) &&
                                 throws(
                                   [&]() { s.cell->computeDistributedLoad( -1, 1, 1, zero, dummy, dummy, 0., 1. ); } ),
                               g.name + ": unknown load types must throw" );
    }

    // consistent tangent, from the undeformed state and after an accepted finite increment
    for ( bool afterStep : { false, true } ) {
      Setup< nDim > s( g );
      if ( afterStep )
        s.accept( increment( s.nDof, nU, 0.03, 3e-3 ) ); // finite deformation and damage history
      const Eigen::VectorXd dQ = afterStep ? increment( s.nDof, nU, 0.01, 1e-3 ) : increment( s.nDof, nU, 0.03, 3e-3 );
      const auto [P, K]        = s.trial( dQ );
      const auto   numK = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ, nU );
      const double err  = ( K - numK ).norm() / numK.norm();
      throwExceptionOnFailure( err < 1e-7,
                               MakeString() << g.name << ( afterStep ? ", second" : ", first" )
                                            << " step: tangent inconsistent, relative error " << err );
      if ( afterStep ) {
        const int    nN    = s.nDof - nU;
        const double errUN = ( K.block( 0, nU, nU, nN ) - numK.block( 0, nU, nU, nN ) ).norm() /
                             numK.block( 0, nU, nU, nN ).norm();
        throwExceptionOnFailure( errUN < 1e-6,
                                 MakeString() << g.name << ": dr_U/dN inconsistent, relative error " << errUN );
      }
    }

    // a finite rigid rotation (about z) does not load the cell
    {
      Setup< nDim >   s( g );
      const double    phi = 0.5;
      Eigen::Matrix3d R   = Eigen::Matrix3d::Identity();
      R.topLeftCorner< 2, 2 >() << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );
      Eigen::VectorXd dQ = Eigen::VectorXd::Zero( s.nDof );
      for ( int A = 0; A < nNodes; A++ ) {
        Eigen::Vector3d X              = Eigen::Vector3d::Zero();
        X.head< nDim >()               = nodeAt( A );
        dQ.segment< nDim >( nDim * A ) = ( ( R - Eigen::Matrix3d::Identity() ) * X ).template head< nDim >();
      }
      const auto P = s.trial( dQ ).first;
      throwExceptionOnFailure( P.norm() < 1e-10,
                               MakeString()
                                 << g.name << ": rigid rotation must not load the cell, |P| = " << P.norm() );
    }

    // body force and lumped mass add up to the totals
    {
      Setup< nDim >   s( g );
      const double    V    = s.mpVolume * s.mps.size();
      const Vec       b    = Vec::LinSpaced( 0.3, -1.1 );
      Eigen::VectorXd fExt = Eigen::VectorXd::Zero( s.nDof );
      Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( s.nDof, s.nDof );
      s.cell->computeBodyLoad( s.cell->getSupportedBodyLoadTypes().at( "BODYFORCE" ),
                               b.data(),
                               fExt.data(),
                               K.data(),
                               0.,
                               1. );
      for ( int i = 0; i < nDim; i++ ) {
        double sum = 0;
        for ( int A = 0; A < nNodes; A++ )
          sum += fExt[nDim * A + i];
        throwExceptionOnFailure( checkIfEqual( sum, -b[i] * V, 1e-12 ), g.name + ": body force must add up to -b V" );
      }

      s.trial( increment( s.nDof, nU, 0.0 ) ); // the points learn their density on the first computation
      Eigen::VectorXd I = Eigen::VectorXd::Zero( s.nDof );
      s.cell->computeLumpedInertia( I.data() );
      throwExceptionOnFailure( checkIfEqual( I.head( nU ).sum(), nDim * matProps[iRho] * V, 1e-12 ),
                               g.name + ": lumped mass must add up to nDim rho V" );
      throwExceptionOnFailure( I.tail( s.nDof - nU ).norm() == 0.0, g.name + ": no inertia on the nonlocal field" );
    }

    // pressure on a material point: consistent follower load tangent
    {
      Setup< nDim >         s( g );
      const Eigen::VectorXd dQ       = increment( s.nDof, nU, 0.03 );
      const Vec             p        = -2.0 * Vec::Unit( nDim - 1 );
      const int             pressure = s.cell->getSupportedDistributedLoadTypes().at( "PRESSURE" );

      auto load = [&]( const Eigen::VectorXd& q ) {
        s.trial( q ); // put the material points into the deformed state
        for ( auto& mp : s.mps )
          mp->prepareYourself( 1.0, 1.0 );
        s.cell->interpolateFieldsToMaterialPoints( q.data() );
        Eigen::VectorXd fExt = Eigen::VectorXd::Zero( s.nDof );
        Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( s.nDof, s.nDof );
        s.cell->computeDistributedLoad( pressure, 1, 3, p.data(), fExt.data(), K.data(), 1.0, 1.0 );
        return std::make_pair( fExt, K );
      };

      const auto [f0, K0] = load( dQ );
      const auto numK     = numericalTangent( [&]( const Eigen::VectorXd& q ) { return load( q ).first; }, dQ, nU );

      throwExceptionOnFailure( f0.norm() > 0, g.name + ": pressure must load the cell" );
      throwExceptionOnFailure( ( K0 - numK ).norm() < 1e-7 * numK.norm(),
                               g.name + ": pressure load tangent inconsistent" );
    }
  }

} // namespace

void testUnknownNamesAreRejected()
{
  const double x[2] = { 0, 0 };
  throwExceptionOnFailure( throws( [&]() {
                             std::unique_ptr< MarmotCell >(
                               MarmotLibrary::MarmotCellFactory::createCell( "NoSuchCell", 1, x, 2 ) );
                           } ),
                           "unknown cell" );
  throwExceptionOnFailure( throws( [&]() {
                             std::unique_ptr< MarmotCell >(
                               MarmotLibrary::MarmotCellFactory::createBSplineCell( "NoSuchCell", 1, x, 2, x, 2 ) );
                           } ),
                           "unknown B-spline cell" );
  throwExceptionOnFailure( throws( [&]() {
                             std::unique_ptr< MarmotMaterialPoint >(
                               MarmotLibrary::MarmotMaterialPointFactory::createMaterialPoint( "NoSuchPoint",
                                                                                               1,
                                                                                               x,
                                                                                               2,
                                                                                               1.0 ) );
                           } ),
                           "unknown material point" );
}

int main()
{
  std::vector< std::function< void() > > testFunctions = { testUnknownNamesAreRejected,
                                                           []() { checkCell< 2 >( lagrangian< 2 >() ); },
                                                           []() { checkCell< 3 >( lagrangian< 3 >() ); } };
  for ( int p = 1; p <= 3; p++ ) {
    testFunctions.push_back( [p]() { checkCell< 2 >( bSpline< 2 >( p ) ); } );
    testFunctions.push_back( [p]() { checkCell< 3 >( bSpline< 3 >( p ) ); } );
  }
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
