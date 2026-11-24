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
#include "Marmot/GradientEnhancedCompressibleNeoHookeDamage.h"
#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotParticleTesting.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Discretizations::Meshfree;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };

  // Both registered materials have dL/dN = 0. This one adds beta * N to the local driving force, so that the
  // dL/dN block of the nonlocal balance is tested as well.
  class LocalFieldDependingOnTheNonlocalField : public Materials::GradientEnhancedCompressibleNeoHookeDamage {
  public:
    static constexpr double beta = 0.3;

    using GradientEnhancedCompressibleNeoHookeDamage::computeStress;
    using GradientEnhancedCompressibleNeoHookeDamage::GradientEnhancedCompressibleNeoHookeDamage;

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override
    {
      GradientEnhancedCompressibleNeoHookeDamage::computeStress( response, tangents, deformation, timeIncrement );
      response.L( 0 ) += beta * deformation.N( 0 );
      tangents.dL_dN( 0, 0 ) += beta;
    }
  };

  const std::string testMaterialName = "TESTLOCALFIELDDEPENDINGONTHENONLOCALFIELD";

  const bool testMaterialRegistered = Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<
    1 >::registerMaterial< LocalFieldDependingOnTheNonlocalField >( testMaterialName );

  // A nonlocal field with a viscosity and a micro-inertia, for the explicit scheme; eta^2/4 = m would be the best m,
  // any positive pair will do here.
  class ViscousNonlocalField : public Materials::GradientEnhancedCompressibleNeoHookeDamage {
  public:
    static constexpr double eta = 0.7;
    static constexpr double m   = 0.05;

    using GradientEnhancedCompressibleNeoHookeDamage::GradientEnhancedCompressibleNeoHookeDamage;

    NonlocalVector getNonlocalViscosity( const double* ) const override { return NonlocalVector::Constant( eta ); }

    NonlocalVector getNonlocalMicroInertia( const double* ) const override { return NonlocalVector::Constant( m ); }
  };

  const std::string viscousMaterialName = "TESTVISCOUSNONLOCALFIELD";

  const bool viscousMaterialRegistered = Marmot::Registry::MaterialGradientEnhancedFiniteStrainFactory<
    1 >::registerMaterial< ViscousNonlocalField >( viscousMaterialName );

  // GM7/COROTATIONAL: {E, nu, k1, k2, k3, k4, k5, c, m, c3, fc, transitionWidth, transitionExponent, stressUnit,
  // density, eta, micro-inertia = eta^2/4}
  const std::vector< double > gm7Props =
    { 2.5e4, 0.18, 1.5e-4, 110, 30, 100, 1e-4, 25.0, 1.0, 4.0, 36.0, 4.0, 0.4, 1.0, 2.4e-9, 3e-4, 2.25e-8 };

  // a grid of 5^nDim 2nd order B-spline kernels with spacing 1 around the particle
  template < int nDim >
  struct KernelGrid {
    std::vector< std::array< double, 3 > >                                                                  centers;
    std::vector< std::unique_ptr< Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed > > kernels;
    std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* >                                 pointers;

    KernelGrid()
    {
      const int nz = nDim == 3 ? 5 : 1;
      for ( int k = 0; k < nz; k++ )
        for ( int j = 0; j < 5; j++ )
          for ( int i = 0; i < 5; i++ )
            centers.push_back( { double( i ), double( j ), double( k ) } );
      for ( auto& c : centers ) {
        kernels.emplace_back(
          std::make_unique< Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed >( c.data(),
                                                                                                     nDim,
                                                                                                     1.6 ) );
        pointers.push_back( kernels.back().get() );
      }
    }
  };

  template < int nDim, int nNonlocal = 1 >
  struct Setup {
    static constexpr int nBlock = nDim + nNonlocal; // interleaved particle dof vector [u_1 N_1 u_2 N_2 ...]

    KernelGrid< nDim >                                                grid;
    Marmot::Discretizations::Meshfree::ReproducingKernelApproximation approximation{ nDim, 1 };
    std::unique_ptr< Marmot::Discretizations::Meshfree::Particle >    particle;
    std::vector< double >                                             stateVars;
    std::vector< int >                                                assignedGridIndices;

    Setup( const std::string&           name,
           const std::vector< double >& vertices,
           double                       volume,
           const std::string&           materialName       = "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
           const std::vector< double >& materialProperties = matProps )
    {
      particle.reset(
        Marmot::Registry::ParticleFactory::createParticle( name,
                                                           1,
                                                           vertices.data(),
                                                           static_cast< int >( vertices.size() ),
                                                           volume,
                                                           materialName,
                                                           materialProperties.data(),
                                                           static_cast< int >( materialProperties.size() ),
                                                           approximation ) );
      stateVars.assign( particle->getNumberOfRequiredStateVars(), 0.0 );
      particle->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      particle->initializeYourself();
      assignKernels();
    }

    // the union over all evaluation points, as the host's particle manager does it
    void assignKernels()
    {
      const int             nEval = particle->getNumberOfEvaluationPoints();
      std::vector< double > eval( nEval * nDim );
      particle->getEvaluationCoordinates( eval.data() );

      std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* > assigned;
      assignedGridIndices.clear();
      for ( size_t g = 0; g < grid.pointers.size(); g++ )
        for ( int e = 0; e < nEval; e++ )
          if ( grid.pointers[g]->isInSupport( &eval[e * nDim] ) ) {
            assigned.push_back( grid.pointers[g] );
            assignedGridIndices.push_back( static_cast< int >( g ) );
            break;
          }
      particle->assignMeshfreeKernelFunctions( assigned );
    }

    int nNodes() const { return static_cast< int >( assignedGridIndices.size() ); }
    int nDof() const { return nNodes() * nBlock; }

    // set properties by name, the others to zero
    void setProperties( const std::map< std::string, double >& values )
    {
      const auto            names = particle->getPropertyNames();
      std::vector< double > p( names.size(), 0.0 );
      for ( size_t i = 0; i < names.size(); i++ )
        if ( values.count( names[i] ) )
          p[i] = values.at( names[i] );
      particle->setProperties( p.data(), static_cast< int >( p.size() ) );
    }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > trial( const Eigen::VectorXd& dQ, double dT = 1.0 )
    {
      const auto      backup = stateVars;
      Eigen::VectorXd P      = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K      = Eigen::MatrixXd::Zero( nDof(), nDof() );
      particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), dT, dT );
      stateVars = backup;
      return { P, K };
    }

    void accept( const Eigen::VectorXd& dQ, double dT = 1.0 )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof(), nDof() );
      particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), dT, dT );
      particle->acceptStateAndPosition();
      assignKernels();
    }

    // a uniform translation t of all nodes, without a nonlocal field
    Eigen::VectorXd translation( const Eigen::Matrix< double, nDim, 1 >& t ) const
    {
      Eigen::VectorXd dQ = Eigen::VectorXd::Zero( nDof() );
      for ( int A = 0; A < nNodes(); A++ )
        dQ.template segment< nDim >( A * nBlock ) = t;
      return dQ;
    }
  };

  template < int nDim, int nNonlocal = 1 >
  Eigen::VectorXd increment( int nNodes, double scaleU, double N0 )
  {
    constexpr int   nBlock = nDim + nNonlocal;
    Eigen::VectorXd dQ( nNodes * nBlock );
    for ( int A = 0; A < nNodes; A++ ) {
      for ( int i = 0; i < nDim; i++ )
        dQ[A * nBlock + i] = scaleU * std::sin( 1.7 * ( A * nDim + i + 1 ) );
      for ( int K = 0; K < nNonlocal; K++ )
        dQ[A * nBlock + nDim + K] = N0 * ( 1. + 0.3 * std::cos( 0.9 * ( A + 1 ) + K ) );
    }
    return dQ;
  }

  // A finite-difference check of a tangent K = df/dQ at Q0: the central differences f'(Q0) D along the columns of D,
  // compared with K D. Up to 64 dofs, D is the identity, i.e. the full finite-difference Jacobian with its column-wise
  // diagnostics; beyond, three fixed random unit directions per field, which exercise the same code with O(1) instead
  // of O(nDof) residual evaluations, so that the tests of the higher-order cells and of the particles with many nodes
  // stay fast also at -O0 (coverage). A wrong tangent fails the directional check with probability one.
  struct NumericalTangent {
    Eigen::MatrixXd     D;        ///< the directions, one per column
    Eigen::MatrixXd     KD;       ///< f'(Q0) D by central differences
    std::vector< bool > nonlocal; ///< direction j acts on the nonlocal dofs only

    /// relative error of K against the difference quotients
    double error( const Eigen::MatrixXd& K ) const { return ( K * D - KD ).norm() / KD.norm(); }

    /// relative error of the rows selected by @p row, along the directions of one field
    double error( const Eigen::MatrixXd& K, const std::function< bool( int ) >& row, bool nonlocalDirections ) const
    {
      const Eigen::MatrixXd KDa = K * D;
      double                num = 0, den = 0;
      for ( int j = 0; j < D.cols(); j++ )
        if ( nonlocal[j] == nonlocalDirections )
          for ( int a = 0; a < KD.rows(); a++ )
            if ( row( a ) ) {
              num += std::pow( KDa( a, j ) - KD( a, j ), 2 );
              den += std::pow( KD( a, j ), 2 );
            }
      return std::sqrt( num / std::max( den, 1e-300 ) );
    }
  };

  NumericalTangent numericalTangentAlong( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                          const Eigen::VectorXd&                                            Q0,
                                          const std::function< bool( int ) >&                               isNonlocal,
                                          double                                                            hU,
                                          double                                                            hN )
  {
    const int        n = static_cast< int >( Q0.size() );
    NumericalTangent t;
    if ( n <= 64 ) {
      t.D = Eigen::MatrixXd::Identity( n, n );
      for ( int j = 0; j < n; j++ )
        t.nonlocal.push_back( isNonlocal( j ) );
    }
    else {
      std::mt19937                             rng( 1 );
      std::uniform_real_distribution< double > uniform( -1.0, 1.0 );
      std::vector< Eigen::VectorXd >           directions;
      for ( bool field : { false, true } ) {
        bool any = false;
        for ( int j = 0; j < n; j++ )
          any = any || isNonlocal( j ) == field;
        for ( int k = 0; any && k < 3; k++ ) {
          Eigen::VectorXd d = Eigen::VectorXd::Zero( n );
          for ( int j = 0; j < n; j++ )
            if ( isNonlocal( j ) == field )
              d[j] = uniform( rng );
          directions.push_back( d.normalized() );
          t.nonlocal.push_back( field );
        }
      }
      t.D.resize( n, directions.size() );
      for ( size_t k = 0; k < directions.size(); k++ )
        t.D.col( k ) = directions[k];
    }
    t.KD.resize( f( Q0 ).size(), t.D.cols() );
    for ( int j = 0; j < t.D.cols(); j++ ) {
      const double h = t.nonlocal[j] ? hN : hU;
      t.KD.col( j )  = ( f( Q0 + h * t.D.col( j ) ) - f( Q0 - h * t.D.col( j ) ) ) / ( 2 * h );
    }
    return t;
  }

  // interleaved layout [u_1 N_1 u_2 N_2 ...]; the nonlocal field is small, hence its smaller step
  template < int nDim >
  NumericalTangent numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                     const Eigen::VectorXd&                                            Q0 )
  {
    return numericalTangentAlong(
      f, Q0, []( int j ) { return ( j % ( nDim + 1 ) ) == nDim; }, 1e-7, 1e-9 );
  }

  // relative errors of the U-U, U-N, N-U, N-N blocks in the interleaved layout
  template < int nDim >
  std::array< double, 4 > blockErrors( const Eigen::MatrixXd& K, const NumericalTangent& numK )
  {
    const auto uRow = []( int a ) { return ( a % ( nDim + 1 ) ) != nDim; };
    const auto nRow = []( int a ) { return ( a % ( nDim + 1 ) ) == nDim; };
    return { numK.error( K, uRow, false ),
             numK.error( K, uRow, true ),
             numK.error( K, nRow, false ),
             numK.error( K, nRow, true ) };
  }

  // the particle geometries: a point, a quad and a hexahedron in the interior of the kernel grid
  std::vector< double > vertices( const std::string& shape )
  {
    if ( shape == "Point" )
      return { 1.9, 2.2 };
    if ( shape == "Point3D" )
      return { 1.9, 2.2, 2.0 };
    if ( shape == "Quad" )
      return { 1.6, 1.9, 2.2, 1.9, 2.2, 2.5, 1.6, 2.5 };
    // Hexa
    std::vector< double > v;
    for ( double z : { 1.7, 2.3 } )
      for ( auto [x, y] :
            std::vector< std::pair< double, double > >{ { 1.6, 1.9 }, { 2.2, 1.9 }, { 2.2, 2.5 }, { 1.6, 2.5 } } ) {
        v.push_back( x );
        v.push_back( y );
        v.push_back( z );
      }
    return v;
  }

  double volume( const std::string& shape )
  {
    return shape == "Hexa" || shape == "Point3D" ? 0.216 : 0.36;
  }

  int numberOfFaces( const std::string& shape )
  {
    return shape == "Quad" ? 4 : shape == "Hexa" ? 6 : 0;
  }

  // fill freed heap chunks of many sizes with garbage, so that members which are read before they are set show up as
  // garbage instead of as zeros or as stale values of a previous, identical object
  void dirtyHeap()
  {
    for ( int n = 1; n <= 1024; n += 1 + n / 8 )
      std::vector< double >( n, 1.2345e300 );
  }

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

  // exact: the full tangent; approximate (NSNI): the stabilization's tangent omits the derivative of
  // dTau/dDeltaF itself (d2tau/dF2 and the damage dependence), so only the N rows are exact there
  enum class Tangent { Exact, NSNIApproximate };

  Tangent tangentMode( const std::string& name )
  {
    return name.find( "NSNI" ) != std::string::npos ? Tangent::NSNIApproximate : Tangent::Exact;
  }

  template < int nDim >
  void checkTangent( const std::string& what, const Eigen::MatrixXd& K, const NumericalTangent& numK, Tangent mode )
  {
    const auto be = blockErrors< nDim >( K, numK );
    if ( mode == Tangent::Exact ) {
      const double err = numK.error( K );
      throwExceptionOnFailure( err < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << what << ": tangent inconsistent, relative error " << err );
    }
    else {
      throwExceptionOnFailure( be[2] < 5e-5 && be[3] < 1e-7,
                               Marmot::Utilities::MakeString()
                                 << what << ": nonlocal rows inconsistent, " << be[2] << ", " << be[3] );
      throwExceptionOnFailure( be[0] < 5e-3,
                               Marmot::Utilities::MakeString()
                                 << what << ": U-U block beyond the NSNI approximation, " << be[0] );
    }
  }

  template < int nDim >
  void checkParticle( const std::string& name, const std::string& shape )
  {
    const auto    v    = vertices( shape );
    const double  V    = volume( shape );
    const auto    mode = tangentMode( name );
    constexpr int nB   = nDim + 1;

    // first step, from the undeformed state
    {
      Setup< nDim >         s( name, v, V );
      const Eigen::VectorXd dQ = increment< nDim >( s.nNodes(), 0.02, 3e-3 );
      const auto [P, K]        = s.trial( dQ );
      const auto numK = numericalTangent< nDim >( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
      checkTangent< nDim >( name + ", first step", K, numK, mode );

      // translation invariance: the displacement residuals of every direction sum to zero
      for ( int i = 0; i < nDim; i++ ) {
        double sum = 0;
        for ( int A = 0; A < s.nNodes(); A++ )
          sum += P[A * nB + i];
        throwExceptionOnFailure( std::abs( sum ) < 1e-10 * P.norm(), name + ": residuals are not self-equilibrated" );
      }
    }

    // second step, from an accepted finite deformation with damage history (and an updated smoothing domain)
    {
      Setup< nDim > s( name, v, V );
      s.accept( increment< nDim >( s.nNodes(), 0.02, 3e-3 ) );
      const Eigen::VectorXd dQ = increment< nDim >( s.nNodes(), 0.01, 1e-3 );
      const auto [P, K]        = s.trial( dQ );
      const auto numK = numericalTangent< nDim >( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
      checkTangent< nDim >( name + ", second step", K, numK, mode );
    }

    // a finite rigid rotation (about z) of the kernel nodes does not load the particle
    {
      Setup< nDim >   s( name, v, V );
      const double    phi = 0.5;
      Eigen::Matrix3d R   = Eigen::Matrix3d::Identity();
      R.topLeftCorner< 2, 2 >() << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );
      Eigen::VectorXd dQ = Eigen::VectorXd::Zero( s.nDof() );
      for ( int A = 0; A < s.nNodes(); A++ ) {
        const Eigen::Vector3d X( s.grid.centers[s.assignedGridIndices[A]].data() );
        dQ.segment< nDim >( A * nB ) = ( ( R - Eigen::Matrix3d::Identity() ) * X ).template head< nDim >();
      }
      const auto P = s.trial( dQ ).first;
      throwExceptionOnFailure( P.norm() < 1e-8,
                               Marmot::Utilities::MakeString()
                                 << name << ": rigid rotation loads the particle, |P| = " << P.norm() );
    }

    // dynamics: with Newmark integration, the inertia enters residual and tangent consistently
    {
      Setup< nDim > s( name, v, V );
      s.setProperties( { { "newmark-beta beta", 0.25 }, { "newmark-beta gamma", 0.5 } } );
      const double          dT = 0.01;
      const Eigen::VectorXd dQ = increment< nDim >( s.nNodes(), 0.003, 3e-3 );
      const auto [P, K]        = s.trial( dQ, dT );
      const auto numK = numericalTangent< nDim >( [&]( const Eigen::VectorXd& q ) { return s.trial( q, dT ).first; },
                                                  dQ );
      checkTangent< nDim >( name + ", dynamic", K, numK, mode );
    }

    // body force per unit undeformed volume: by the partition of unity of the test functions, the displacement loads
    // sum to -b V0 (the host's sign convention); the nonlocal rows are not loaded, and a dead load has no tangent
    {
      Setup< nDim > s( name, v, V );
      s.trial( increment< nDim >( s.nNodes(), 0.03, 3e-3 ) );
      const double    b[3] = { 1.5, -2.0, 0.7 };
      Eigen::VectorXd P    = Eigen::VectorXd::Zero( s.nDof() );
      Eigen::MatrixXd K    = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
      s.particle
        ->computeBodyLoad( s.particle->getSupportedBodyLoadTypes().at( "BODYFORCE" ), b, P.data(), K.data(), 1.0, 1.0 );
      const double V0 = s.particle->getVolumeUndeformed();
      for ( int i = 0; i < nDim; i++ ) {
        double sum = 0;
        for ( int A = 0; A < s.nNodes(); A++ )
          sum += P[A * nB + i];
        throwExceptionOnFailure( std::abs( sum + b[i] * V0 ) < 1e-12 * ( 1. + V0 ),
                                 Marmot::Utilities::MakeString()
                                   << name << ": body force " << sum << " != -b V0 = " << -b[i] * V0 );
      }
      for ( int A = 0; A < s.nNodes(); A++ )
        throwExceptionOnFailure( P[A * nB + nDim] == 0.0, name + ": the body force loads the nonlocal field" );
      throwExceptionOnFailure( K.norm() == 0.0, name + ": a dead body load has no tangent" );
    }

    // a point particle has no faces: no distributed loads
    if ( numberOfFaces( shape ) == 0 ) {
      Setup< nDim > s( name, v, V );
      throwExceptionOnFailure( s.particle->getSupportedDistributedLoadTypes().empty(),
                               name + ": a point particle advertises distributed loads" );
      Eigen::VectorXd P     = Eigen::VectorXd::Zero( s.nDof() );
      Eigen::MatrixXd K     = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
      const double    p     = 1.0;
      bool            threw = false;
      try {
        s.particle->computeDistributedLoad( 0, 1, &p, P.data(), K.data(), 1.0, 1.0 );
      }
      catch ( const std::invalid_argument& ) {
        threw = true;
      }
      throwExceptionOnFailure( threw, name + ": a distributed load on a point particle must throw" );
    }

    // distributed loads (pressure and the correction of the weak form) on a face, with their tangents
    if ( numberOfFaces( shape ) > 0 ) {
      for ( const auto& [loadName, type] : Setup< nDim >( name, v, V ).particle->getSupportedDistributedLoadTypes() ) {
        Setup< nDim >         s( name, v, V );
        const double          p    = 2.0;
        const Eigen::VectorXd dQ   = increment< nDim >( s.nNodes(), 0.03, 3e-3 );
        auto                  load = [&]( const Eigen::VectorXd& q ) {
          s.trial( q ); // the loads read the configuration of the last computation
          Eigen::VectorXd P  = Eigen::VectorXd::Zero( s.nDof() );
          Eigen::MatrixXd K  = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
          Eigen::VectorXd Pi = P;
          Eigen::MatrixXd Ki = K;
          s.particle->computePhysicsKernels( q.data(), Pi.data(), Ki.data(), 1.0, 1.0 );
          s.particle->computeDistributedLoad( type, 2, &p, P.data(), K.data(), 1.0, 1.0 );
          return std::make_pair( P, K );
        };
        const auto backup   = s.stateVars;
        const auto [P0, K0] = load( dQ );
        const auto numK     = numericalTangent< nDim >(
          [&]( const Eigen::VectorXd& q ) {
            const auto r = load( q ).first;
            s.stateVars  = backup;
            return r;
          },
          dQ );
        s.stateVars = backup;
        throwExceptionOnFailure( P0.norm() > 0, name + ": " + loadName + " must load the particle" );
        const double err = std::min( numK.error( K0 ), numK.error( -K0 ) );
        throwExceptionOnFailure( err < 1e-6,
                                 Marmot::Utilities::MakeString()
                                   << name << ": " << loadName << " load tangent inconsistent, " << err );
        for ( int A = 0; A < s.nNodes(); A++ )
          throwExceptionOnFailure( P0[A * nB + nDim] == 0,
                                   name + ": " + loadName + " must not load the nonlocal field" );

        // as for the elements, a positive pressure pushes on the face: in the undeformed configuration, the residual
        // of the load (P = -f_ext, summed over the nodes by the partition of unity) is p |face| n, with the outward
        // normal n of the axis-aligned face pointing from the particle center to the face center
        if ( loadName == "PRESSURE" ) {
          using Vec = Eigen::Matrix< double, nDim, 1 >;
          Vec center, faceCenter, PSum = Vec::Zero();
          s.particle->getCenterCoordinates( center.data() );
          s.particle->getFaceCoordinates( 2, faceCenter.data() );
          const Vec    n           = ( faceCenter - center ).normalized();
          const double faceArea    = std::pow( 0.6, nDim - 1 );
          const auto   PUndeformed = load( Eigen::VectorXd::Zero( s.nDof() ) ).first;
          s.stateVars              = backup;
          for ( int A = 0; A < s.nNodes(); A++ )
            PSum += PUndeformed.template segment< nDim >( A * nB );
          throwExceptionOnFailure( ( PSum - p * faceArea * n ).norm() < 1e-10 * p * faceArea,
                                   Marmot::Utilities::MakeString()
                                     << name << ": a positive pressure must push on the face, P = " << PSum.transpose()
                                     << " != p |face| n = " << ( p * faceArea * n ).transpose() );
        }
      }
    }
  }

  // the weak-form correction on all faces balances the displacement residual of a homogeneous deformation
  // (divergence theorem of the smoothed gradient); after a step, only the full SQCNI smooths over the deformed
  // particle itself, which makes the balance exact
  // the tangent including dL/dN, with a material whose local driving force depends on the nonlocal field
  template < int nDim >
  void checkTangentWithdL_dN( const std::string& name, const std::string& shape )
  {
    Setup< nDim >         s( name, vertices( shape ), volume( shape ), testMaterialName );
    const Eigen::VectorXd dQ = increment< nDim >( s.nNodes(), 0.02, 3e-3 );
    const auto [P, K]        = s.trial( dQ );
    const auto numK = numericalTangent< nDim >( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
    checkTangent< nDim >( name + ", dL/dN != 0", K, numK, tangentMode( name ) );
  }

  template < int nDim >
  void checkWeakFormCorrection( const std::string& name, const std::string& shape )
  {
    using Vec          = Eigen::Matrix< double, nDim, 1 >;
    constexpr int nB   = nDim + 1;
    const auto    v    = vertices( shape );
    const bool    full = name.find( "SQCNI/" ) != std::string::npos || name.find( "SQCNIxNSNI/" ) != std::string::npos;

    for ( bool afterStep : { false, true } ) {
      if ( afterStep && !full )
        continue;
      Setup< nDim >                       s( name, v, volume( shape ) );
      Eigen::Matrix< double, nDim, nDim > G = Eigen::Matrix< double, nDim, nDim >::Zero();
      G.diagonal().setConstant( 0.04 );
      G( 0, 1 )        = 0.03;
      auto homogeneous = [&]( const Eigen::Matrix< double, nDim, nDim >& H ) {
        Eigen::VectorXd dQ = Eigen::VectorXd::Zero( s.nDof() );
        for ( int A = 0; A < s.nNodes(); A++ )
          dQ.template segment< nDim >( A * nB ) = H * Eigen::Map< const Vec >(
                                                        s.grid.centers[s.assignedGridIndices[A]].data() );
        return dQ;
      };
      if ( afterStep )
        s.accept( homogeneous( G ) );
      const Eigen::VectorXd dQ = homogeneous( 0.5 * G.transpose() );

      // the loads read the configuration of the last computation
      Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() ), Pc = P;
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() ), Kc = K;
      s.particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
      const int cwf = s.particle->getSupportedDistributedLoadTypes().at( "CWFCORRECTION" );
      for ( int f = 1; f <= numberOfFaces( shape ); f++ )
        s.particle->computeDistributedLoad( cwf, f, nullptr, Pc.data(), Kc.data(), 1.0, 1.0 );

      double res = 0, ref = 0;
      for ( int A = 0; A < s.nNodes(); A++ ) {
        res += ( P + Pc ).template segment< nDim >( A * nB ).squaredNorm();
        ref += P.template segment< nDim >( A * nB ).squaredNorm();
      }
      throwExceptionOnFailure( std::sqrt( res ) < 1e-10 * std::sqrt( ref ),
                               Marmot::Utilities::MakeString()
                                 << name << ( afterStep ? ", second step" : ", first step" )
                                 << ": weak-form correction does not balance a homogeneous state, "
                                 << std::sqrt( res / ref ) );

      // a component mask (load 1 = x only) keeps exactly those rows of the full correction, residual and tangent
      {
        Eigen::VectorXd Pm    = Eigen::VectorXd::Zero( s.nDof() );
        Eigen::MatrixXd Km    = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
        const double    xOnly = 1.0;
        for ( int f = 1; f <= numberOfFaces( shape ); f++ )
          s.particle->computeDistributedLoad( cwf, f, &xOnly, Pm.data(), Km.data(), 1.0, 1.0 );
        double err = 0;
        for ( int r = 0; r < s.nDof(); r++ ) {
          const double kept = r % nB == 0 ? 1.0 : 0.0;
          err += std::abs( Pm( r ) - kept * Pc( r ) ) + ( Km.row( r ) - kept * Kc.row( r ) ).cwiseAbs().sum();
        }
        throwExceptionOnFailure( err < 1e-12 * ( 1 + Pc.norm() + Kc.norm() ),
                                 Marmot::Utilities::MakeString()
                                   << name << ": masked weak-form correction differs from the selected rows, " << err );
      }
    }
  }

  // the parts of the particle interface beyond the mechanics: properties, geometry, state, VCI
  template < int nDim >
  void checkInterface( const std::string& name, const std::string& shape )
  {
    using Vec        = Eigen::Matrix< double, nDim, 1 >;
    const auto    v  = vertices( shape );
    const double  V  = volume( shape );
    constexpr int nB = nDim + 1;

    // properties, fields, dimension
    {
      Setup< nDim > s( name, v, V );
      const auto    names = s.particle->getPropertyNames();
      throwExceptionOnFailure( std::find( names.begin(), names.end(), "VCI order" ) != names.end(),
                               name + ": VCI order is a property" );
      std::vector< double > tooFew( names.size() - 1, 0.0 );
      throwExceptionOnFailure( throws( [&]() {
                                 s.particle->setProperties( tooFew.data(), static_cast< int >( tooFew.size() ) );
                               } ),
                               name + ": a wrong number of properties must throw" );
      const double one = 1.0;
      throwExceptionOnFailure( throws( [&]() { s.particle->setProperty( "no such property", &one ); } ),
                               name + ": an unknown property must throw" );

      s.setProperties( { { "VCI order", 1.0 } } );
      throwExceptionOnFailure( s.particle->vci_getNumberOfConstraints() == nDim + 1,
                               name + ": number of linear VCI constraints" );

      throwExceptionOnFailure( s.particle->getFields() ==
                                 std::vector< std::string >{ "displacement", "nonlocal damage" },
                               name + ": fields" );
      throwExceptionOnFailure( s.particle->getNBaseDof() == nB && s.particle->getDimension() == nDim,
                               name + ": dimension" );
      throwExceptionOnFailure( !s.particle->getParticleShape().empty(), name + ": shape" );
      throwExceptionOnFailure( throws( [&]() { s.particle->setInitialCondition( "no such condition", &one ); } ),
                               name + ": an unknown initial condition must throw" );
      throwExceptionOnFailure( s.particle->getSupportedBodyLoadTypes().count( "BODYFORCE" ) == 1,
                               name + ": body force type" );
    }

    // geometry of the undeformed particle
    {
      Setup< nDim > s( name, v, V );
      const int     nV = s.particle->getNumberOfVertices();
      throwExceptionOnFailure( nV * nDim == int( v.size() ), name + ": number of vertices" );

      Eigen::VectorXd vert( nV * nDim ), vis( nV * nDim );
      s.particle->getVertexCoordinates( vert.data() );
      s.particle->getVisualizationVertexCoordinates( vis.data() );
      const Eigen::Map< const Eigen::VectorXd > v0( v.data(), v.size() );
      throwExceptionOnFailure( ( vert - v0 ).norm() < 1e-14 && ( vis - v0 ).norm() < 1e-14,
                               name + ": undeformed vertices" );

      Vec center, mean = Vec::Zero();
      s.particle->getCenterCoordinates( center.data() );
      for ( int k = 0; k < nV; k++ )
        mean += v0.template segment< nDim >( k * nDim ) / nV;
      throwExceptionOnFailure( ( center - mean ).norm() < 1e-14, name + ": center" );

      if ( numberOfFaces( shape ) == 0 )
        throwExceptionOnFailure( throws( [&]() { s.particle->getFaceCoordinates( 1, center.data() ); } ),
                                 name + ": a point has no faces" );
      else {
        Vec faceMean = Vec::Zero();
        for ( int f = 1; f <= numberOfFaces( shape ); f++ ) {
          Vec c;
          s.particle->getFaceCoordinates( f, c.data() );
          faceMean += c / numberOfFaces( shape );
        }
        throwExceptionOnFailure( ( faceMean - mean ).norm() < 1e-14, name + ": face centers" );
      }

      Eigen::VectorXd N( s.nNodes() );
      s.particle->getInterpolationVector( N.data(), center.data() );
      throwExceptionOnFailure( std::abs( N.sum() - 1.0 ) < 1e-12, name + ": interpolation is a partition of unity" );
    }

    // a translation is carried by the material point, the center and the vertices, and loads nothing
    {
      Setup< nDim > s( name, v, V );
      const Vec     t = Vec::LinSpaced( 0.1, 0.1 * nDim );
      s.accept( s.translation( t ) );

      const auto u = s.particle->getStateView( "displacement", 0 );
      throwExceptionOnFailure( ( Eigen::Map< const Vec >( u.stateLocation ) - t ).norm() < 1e-12,
                               name + ": displacement of a translation" );
      if ( numberOfFaces( shape ) > 0 ) {
        const auto view = s.particle->getStateView( "vertex displacements", 0 );
        throwExceptionOnFailure( view.stateSize == int( v.size() ), name + ": size of the vertex displacements" );
        for ( size_t k = 0; k < v.size() / nDim; k++ )
          throwExceptionOnFailure( ( Eigen::Map< const Vec >( view.stateLocation + k * nDim ) - t ).norm() < 1e-12,
                                   name + ": vertex displacements of a translation" );
      }

      Vec translated, mean = Vec::Zero();
      s.particle->getCenterCoordinates( translated.data() );
      for ( size_t k = 0; k < v.size() / nDim; k++ )
        mean += Eigen::Map< const Vec >( &v[k * nDim] ) / ( v.size() / nDim );
      throwExceptionOnFailure( ( translated - mean - t ).norm() < 1e-12, name + ": translated center" );
    }

    // geostatic initial stress: the material point accepts it
    {
      Setup< nDim > s( name, v, V );
      const double  geostatic[5] = { -1.0, 1.0, -1.0, 1.0, 0.5 };
      s.particle->setInitialCondition( "geostaticstress", geostatic );
    }

    // VCI of order 0, before the first increment: the correction eta = M^-1 R makes the corrected test functions
    // satisfy R = int_dOmega T n - int_Omega grad T = 0 for every node whose correction is active (M > 0)
    if ( numberOfFaces( shape ) > 0 ) {
      dirtyHeap(); // so that an unevaluated basis or volume shows up as garbage
      Setup< nDim > s( name, v, V );
      const int     nR = s.nNodes() * nDim;
      throwExceptionOnFailure( s.particle->vci_getNumberOfConstraints() == 1, name + ": one constraint of order 0" );

      auto residual = [&]() {
        Eigen::VectorXd R = Eigen::VectorXd::Zero( nR ), vol = Eigen::VectorXd::Zero( nR );
        for ( int f = 1; f <= numberOfFaces( shape ); f++ )
          s.particle->vci_compute_Test_P_BoundaryIntegral( R.data(), nullptr, f );
        s.particle->vci_compute_TestGradient_P_Integral( vol.data() );
        s.particle->vci_compute_Test_PGradient_Integral( vol.data() );
        return Eigen::VectorXd( R - vol );
      };

      Eigen::VectorXd M = Eigen::VectorXd::Zero( s.nNodes() );
      s.particle->vci_compute_MMatrix( M.data() );
      throwExceptionOnFailure( std::abs( M.maxCoeff() - V ) < 1e-12,
                               Marmot::Utilities::MakeString()
                                 << name << ": M of order 0 is the volume, " << M.maxCoeff() << " vs " << V );

      const Eigen::VectorXd R0 = residual();
      Eigen::VectorXd       eta( nR );
      for ( int A = 0; A < s.nNodes(); A++ )
        eta.template segment< nDim >( A * nDim ) = M[A] > 0
                                                     ? Eigen::VectorXd( R0.template segment< nDim >( A * nDim ) / M[A] )
                                                     : Eigen::VectorXd::Zero( nDim );
      s.particle->vci_assignTestFunctionCorrectionTerms( eta.data() );
      const Eigen::VectorXd R1 = residual();
      for ( int A = 0; A < s.nNodes(); A++ )
        if ( M[A] > 0 )
          throwExceptionOnFailure( R1.template segment< nDim >( A * nDim ).norm() < 1e-12,
                                   Marmot::Utilities::MakeString() << name << ": VCI residual after the correction, "
                                                                   << R1.template segment< nDim >( A * nDim ).norm() );
    }
    else {
      // a point particle contributes T n dA at its center to the boundary integral
      Setup< nDim >   s( name, v, V );
      Eigen::VectorXd R = Eigen::VectorXd::Zero( s.nNodes() * nDim );
      const Vec       n = Vec::LinSpaced( 0.3, 0.6 );
      s.particle->vci_compute_Test_P_BoundaryIntegral( R.data(), n.data(), 0 );
      Vec sum = Vec::Zero();
      for ( int A = 0; A < s.nNodes(); A++ )
        sum += R.template segment< nDim >( A * nDim );
      throwExceptionOnFailure( ( sum - n ).norm() < 1e-12, name + ": point boundary integral" );
    }
  }

  // A point of the particle's own (a vertex, the center) moves with getInterpolationVectorLinearizedAboutCenter():
  // after an accepted increment dQ from a deformed state, x_new - x = sum_B G_B dQ^u_B to round-off, and G sums to one,
  // also for the smoothing domains that do not follow F (the vertices are the geometry, not the smoothing domain). The
  // shape functions at a vertex, getInterpolationVector(), do not move it, the reason a constraint that measures a
  // vertex must distribute its force with G.
  template < int nDim >
  void checkOwnPointsMoveWithTheLinearizedInterpolation( const std::string& name, const std::string& shape )
  {
    using Points          = Eigen::Matrix< double, nDim, Eigen::Dynamic >;
    constexpr int nBlock  = Setup< nDim >::nBlock;
    const bool    isPoint = numberOfFaces( shape ) == 0;

    // a host may ask for the vertices before it initializes the particle (EdelweissMeshfree creates its constraints
    // before it prepares the model): on a zero state block they are the undeformed ones
    {
      const auto                                                        v = vertices( shape );
      Marmot::Discretizations::Meshfree::ReproducingKernelApproximation approximation{ nDim, 1 };
      std::unique_ptr< Marmot::Discretizations::Meshfree::Particle >    particle(
        Marmot::Registry::ParticleFactory::createParticle( name,
                                                           1,
                                                           v.data(),
                                                           static_cast< int >( v.size() ),
                                                           volume( shape ),
                                                           "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                                           matProps.data(),
                                                           static_cast< int >( matProps.size() ),
                                                           approximation ) );
      std::vector< double > stateVars( particle->getNumberOfRequiredStateVars(), 0.0 );
      particle->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      Eigen::VectorXd vert( v.size() );
      particle->getVertexCoordinates( vert.data() );
      throwExceptionOnFailure( ( vert - Eigen::Map< const Eigen::VectorXd >( v.data(), v.size() ) ).norm() < 1e-14,
                               name + ": the vertices before initializeYourself() are the undeformed ones" );
    }

    Setup< nDim > s( name, vertices( shape ), volume( shape ) );
    s.accept( increment< nDim >( s.nNodes(), 0.05, 1e-3 ) ); // a deformed, accepted state

    const int nV = s.particle->getNumberOfVertices();
    const int nN = s.nNodes();
    Points    points( nDim, nV + 1 );
    s.particle->getVertexCoordinates( points.data() );
    s.particle->getCenterCoordinates( points.col( nV ).data() );

    Eigen::MatrixXd G( nN, nV + 1 ), N( nN, nV + 1 );
    for ( int k = 0; k <= nV; k++ ) {
      s.particle->getInterpolationVectorLinearizedAboutCenter( G.col( k ).data(), points.col( k ).data() );
      s.particle->getInterpolationVector( N.col( k ).data(), points.col( k ).data() );
    }

    Eigen::VectorXd dQ( s.nDof() );
    for ( int i = 0; i < s.nDof(); i++ )
      dQ[i] = 0.03 * std::cos( 2.3 * ( i + 1 ) );
    const Eigen::Map< const Points, 0, Eigen::OuterStride<> > dq( dQ.data(), nDim, nN, Eigen::OuterStride<>( nBlock ) );
    s.accept( dQ );

    Points moved( nDim, nV + 1 );
    s.particle->getVertexCoordinates( moved.data() );
    s.particle->getCenterCoordinates( moved.col( nV ).data() );
    const Points motion = moved - points;

    throwExceptionOnFailure( ( motion - dq * G ).norm() < 1e-12 * motion.norm(),
                             name +
                               ": a point of the particle moves with getInterpolationVectorLinearizedAboutCenter()" );
    throwExceptionOnFailure( ( G.colwise().sum().array() - 1.0 ).abs().maxCoeff() < 1e-12,
                             name + ": the linearized interpolation is a partition of unity" );
    if ( !isPoint )
      throwExceptionOnFailure( ( motion.leftCols( nV ) - dq * N.leftCols( nV ) ).norm() > 1e-2 * motion.norm(),
                               name + ": the shape functions at the vertices do not move them" );
  }

} // namespace

void testSQCNIxNSNIFirstIncrementIsDeterministic()
{
  // regression: the intermediate second moments used to be set only on the first accepted increment, so the
  // stabilization of the very first increment was scaled by whatever the heap held. Two particles built after
  // different heap histories must agree bit for bit.
  const std::string name = "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad";

  Setup< 2 > a( name, vertices( "Quad" ), 0.36 );
  const auto Pa = a.trial( increment< 2 >( a.nNodes(), 0.02, 3e-3 ) ).first;

  {
    std::vector< double > garbage( 1 << 16, 1.2345e300 ); // dirty the heap, then free it for reuse
  }
  Setup< 2 > b( name, vertices( "Quad" ), 0.36 );
  const auto Pb = b.trial( increment< 2 >( b.nNodes(), 0.02, 3e-3 ) ).first;

  throwExceptionOnFailure( Pa.allFinite() && ( Pa - Pb ).norm() == 0.0,
                           "first-increment residual depends on uninitialized memory" );
}

namespace {

  // The explicit internal force of the first increment from the virgin state equals the implicit residual without
  // inertia: the nonlocal balance in total form equals the one in increment form when the state it starts from is zero,
  // and both evaluate the same material point.
  template < int nDim, int nNonlocal >
  void checkExplicitAgreesWithImplicit( const std::string&           name,
                                        const std::string&           shape,
                                        const std::string&           materialName,
                                        const std::vector< double >& properties,
                                        double                       scaleU )
  {
    Setup< nDim, nNonlocal > s( name, vertices( shape ), volume( shape ), materialName, properties );
    const Eigen::VectorXd    dQ = increment< nDim, nNonlocal >( s.nNodes(), scaleU, 1e-4 );

    const Eigen::VectorXd PImplicit = s.trial( dQ ).first;

    s.particle->updatePhysicsExplicit( dQ.data(), 1.0, 1.0 );
    Eigen::VectorXd PExplicit = Eigen::VectorXd::Zero( s.nDof() );
    s.particle->computePhysicsKernelsExplicit( PExplicit.data() );

    throwExceptionOnFailure( ( PExplicit - PImplicit ).cwiseAbs().maxCoeff() <= 1e-12 * PImplicit.cwiseAbs().maxCoeff(),
                             Marmot::Utilities::MakeString()
                               << name << ": explicit force differs from the implicit residual by "
                               << ( PExplicit - PImplicit ).cwiseAbs().maxCoeff() << " of "
                               << PImplicit.cwiseAbs().maxCoeff() );
  }

  // Lumped inertia, damping and momentum against hand values: rho0 T_A V0 and m T_A V0, eta T_A V0 on the nonlocal
  // rows only, and the inertia times the rates of a uniform increment (translation and nonlocal field).
  template < int nDim, int nNonlocal >
  void checkLumpedOperators( const std::string&           name,
                             const std::string&           shape,
                             const std::string&           materialName,
                             const std::vector< double >& properties,
                             double                       density,
                             double                       eta,
                             double                       microInertia )
  {
    Setup< nDim, nNonlocal > s( name, vertices( shape ), volume( shape ), materialName, properties );
    constexpr int            nBlock = nDim + nNonlocal;
    const double             V0     = s.particle->getVolumeUndeformed();

    std::array< double, nDim > center;
    s.particle->getCenterCoordinates( center.data() );
    Eigen::VectorXd T( s.nNodes() );
    s.particle->getInterpolationVector( T.data(), center.data() );

    // a uniform increment over dT: v = t / dT, dN/dt = dn / dT, and no gradients
    const double    dT = 1e-3, dn = 2e-6;
    Eigen::VectorXd dQ = Eigen::VectorXd::Zero( s.nDof() );
    for ( int A = 0; A < s.nNodes(); A++ ) {
      for ( int i = 0; i < nDim; i++ )
        dQ( A * nBlock + i ) = 1e-6 * ( i + 1 );
      for ( int K = 0; K < nNonlocal; K++ )
        dQ( A * nBlock + nDim + K ) = dn;
    }
    s.particle->updatePhysicsExplicit( dQ.data(), dT, dT );

    Eigen::VectorXd M = Eigen::VectorXd::Zero( s.nDof() ), C = M, P = M;
    s.particle->computeLumpedInertia( M.data() );
    s.particle->computeLumpedDamping( C.data() );
    s.particle->computeLumpedMomentum( P.data() );

    const auto close = [&]( double a, double b, double scale ) { return std::abs( a - b ) <= 1e-12 * scale; };
    for ( int A = 0; A < s.nNodes(); A++ ) {
      for ( int i = 0; i < nDim; i++ ) {
        const int row = A * nBlock + i;
        throwExceptionOnFailure( close( M( row ), density * T( A ) * V0, density * V0 ), name + ": lumped mass" );
        throwExceptionOnFailure( C( row ) == 0.0, name + ": the displacement rows must have no damping" );
        throwExceptionOnFailure( close( P( row ),
                                        density * T( A ) * V0 * 1e-6 * ( i + 1 ) / dT,
                                        density * V0 * 1e-6 * nDim / dT ),
                                 name + ": lumped momentum, displacement" );
      }
      for ( int K = 0; K < nNonlocal; K++ ) {
        const int row = A * nBlock + nDim + K;
        throwExceptionOnFailure( close( M( row ), microInertia * T( A ) * V0, microInertia * V0 ),
                                 name + ": lumped micro-inertia" );
        throwExceptionOnFailure( close( C( row ), eta * T( A ) * V0, eta * V0 ), name + ": lumped damping" );
        throwExceptionOnFailure( close( P( row ), microInertia * T( A ) * V0 * dn / dT, microInertia * V0 * dn / dT ),
                                 name + ": lumped momentum, nonlocal" );
      }
    }

    // and the totals: the partition of unity distributes exactly the particle's mass and momentum
    double totalMass = 0.0;
    for ( int A = 0; A < s.nNodes(); A++ )
      totalMass += M( A * nBlock );
    throwExceptionOnFailure( close( totalMass, density * V0, density * V0 ), name + ": total mass" );
  }

  // The bound of |K| W of the particle for the critical time step (addBoundOfAbsoluteStiffnessTimes), against
  // the stiffness itself, the tangent of a zero increment, at a state off the reference: never below |K| W, and its
  // largest entry of each column (what the Collatz-Wielandt bound takes) at most maxOverestimate times the exact one
  // (the Cauchy-Schwarz inequality over the nodes overestimates the gradient terms by about 1.2 in 2D, 1.3 in 3D).
  // The weights are positive and varying on the dofs of one field and zero on the other, as the host's field blocks.
  template < int nDim, int nNonlocal >
  void checkBoundOfAbsoluteStiffness( const std::string&           name,
                                      const std::string&           shape,
                                      const std::string&           materialName,
                                      const std::vector< double >& properties,
                                      double                       scaleU,
                                      double                       maxOverestimate )
  {
    constexpr int            nBlock = nDim + nNonlocal;
    Setup< nDim, nNonlocal > s( name, vertices( shape ), volume( shape ), materialName, properties );
    s.accept( increment< nDim, nNonlocal >( s.nNodes(), scaleU, 1e-2 * scaleU ) );

    const int             n = s.nDof();
    const Eigen::MatrixXd K = s.trial( Eigen::VectorXd::Zero( n ) ).second;
    Eigen::MatrixXd       W = Eigen::MatrixXd::Zero( n, 2 );
    for ( int a = 0; a < n; a++ )
      W( a, a % nBlock < nDim ? 0 : 1 ) = 1.0 + 0.5 * std::sin( 0.7 * a );

    Marmot::Testing::checkAbsoluteStiffnessBoundAgainstStiffness( name, *s.particle, K, W, maxOverestimate );

    // the description of the stiffness leaves the particle as it was, also in the middle of an increment
    const Eigen::VectorXd dQ = increment< nDim, nNonlocal >( s.nNodes(), 0.5 * scaleU, 0.5e-2 * scaleU );
    Marmot::Testing::checkDescribingTheStiffnessLeavesTheParticleUnchanged( name, *s.particle, s.stateVars, dQ, true );
  }

  // The explicit force after an accepted increment, where F_n != I and the nonlocal field is no longer zero: the
  // nonlocal rows against r^N_A = (T_A (N - L) + c grad_X T_A . grad_X N) V0 by hand, with grad_X T_A = F_n^T grad_Y
  // T_A and the total grad_X N = grad_Y dN_1 + F_n^T grad_Y dN_2 of the two increments, each from the shape functions
  // of its own configuration. A point particle, whose test functions are the shape functions of the approximation.
  template < int nDim >
  void checkExplicitTotalFormAfterAnAcceptedIncrement( const std::string& name, const std::string& shape )
  {
    Setup< nDim > s( name, vertices( shape ), volume( shape ) );
    constexpr int nBlock = nDim + 1;
    const double  c      = matProps[4] * matProps[4];
    const double  V0     = s.particle->getVolumeUndeformed();

    // the shape functions and their gradients at the particle center, with the kernels currently assigned
    const auto shapeFunctions = [&]() {
      std::array< double, nDim > center;
      s.particle->getCenterCoordinates( center.data() );
      std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* > kernels;
      for ( int g : s.assignedGridIndices )
        kernels.push_back( s.grid.pointers[g] );
      Eigen::VectorXd N( s.nNodes() );
      Eigen::MatrixXd dN_dY( nDim, s.nNodes() );
      s.approximation.computeShapeFunctionsAndGradients( center.data(), kernels, N.data(), dN_dY.data() );
      return std::make_pair( N, dN_dY );
    };
    // grad_Y of the nonlocal increment dQ
    const auto nonlocalGradient = [&]( const Eigen::MatrixXd& dN_dY, const Eigen::VectorXd& dQ ) {
      Eigen::Matrix< double, nDim, 1 > g = Eigen::Matrix< double, nDim, 1 >::Zero();
      for ( int B = 0; B < s.nNodes(); B++ )
        g += dQ( B * nBlock + nDim ) * dN_dY.col( B );
      return g;
    };

    const Eigen::VectorXd            dQ1    = increment< nDim >( s.nNodes(), 0.02, 1e-4 );
    const auto                       dN_dY1 = shapeFunctions().second;
    Eigen::Matrix< double, nDim, 1 > dN_dX  = nonlocalGradient( dN_dY1, dQ1 ); // F_0 = I
    s.particle->updatePhysicsExplicit( dQ1.data(), 1.0, 1.0 );
    s.particle->acceptStateAndPosition();
    s.assignKernels();

    Eigen::Matrix< double, nDim, nDim > Fn;
    const double* F = s.particle->getStateView( "deformation gradient", 0 ).stateLocation; // 3x3, row major
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        Fn( i, j ) = F[3 * i + j];
    throwExceptionOnFailure( ( Fn - Fn.transpose() ).norm() > 1e-4, name + ": the test needs a non-symmetric F_n" );

    const Eigen::VectorXd dQ2 = increment< nDim >( s.nNodes(), 0.01, 5e-5 );
    const auto [N2, dN_dY2]   = shapeFunctions();
    dN_dX += Fn.transpose() * nonlocalGradient( dN_dY2, dQ2 );
    s.particle->updatePhysicsExplicit( dQ2.data(), 2.0, 1.0 );
    Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() );
    s.particle->computePhysicsKernelsExplicit( P.data() );

    const double N     = *s.particle->getStateView( "nonlocal damage", 0 ).stateLocation;
    const double L     = *s.particle->getStateView( "local damage", 0 ).stateLocation;
    const double scale = ( std::abs( N ) + std::abs( L ) + c * dN_dX.norm() ) * V0;
    for ( int A = 0; A < s.nNodes(); A++ ) {
      const double expected = ( N2( A ) * ( N - L ) + c * ( Fn.transpose() * dN_dY2.col( A ) ).dot( dN_dX ) ) * V0;
      throwExceptionOnFailure( std::abs( P( A * nBlock + nDim ) - expected ) <= 1e-10 * scale,
                               Marmot::Utilities::MakeString()
                                 << name << ": nonlocal row " << A << " of the explicit force "
                                 << P( A * nBlock + nDim ) << " differs from the total form " << expected );
    }
  }

  // The lumped momentum of a linear increment dQ_A = G Y_A over dT, i.e. v = G Y / dT: its sum over the nodes is the
  // momentum rho0 V0 G Y_c / dT of the particle, and its first moment sum_A P_A (x) Y_A adds to rho0 V0 v_c (x) Y_c the
  // NSNI stabilization rho0 dv/dY M^Y (with sum_A dT_A/dY (x) Y_A = I, and J_Y = 1 before the first accept), with the
  // second moments M^Y = V0 a^2/12 I of the cube or square of edge a.
  template < int nDim, int nNonlocal >
  void checkLumpedMomentumOfALinearIncrement( const std::string&           name,
                                              const std::string&           shape,
                                              const std::string&           materialName,
                                              const std::vector< double >& properties,
                                              double                       density,
                                              bool                         withNSNIStabilization )
  {
    Setup< nDim, nNonlocal > s( name, vertices( shape ), volume( shape ), materialName, properties );
    constexpr int            nBlock = nDim + nNonlocal;
    const double             V0     = s.particle->getVolumeUndeformed();
    const double             dT     = 1e-3;

    Eigen::Matrix< double, nDim, nDim > G;
    for ( int i = 0; i < nDim; i++ )
      for ( int j = 0; j < nDim; j++ )
        G( i, j ) = 1e-6 * std::sin( 1.1 * ( nDim * i + j + 1 ) );

    Eigen::Matrix< double, nDim, Eigen::Dynamic > Y( nDim, s.nNodes() );
    Eigen::VectorXd                               dQ = Eigen::VectorXd::Zero( s.nDof() );
    for ( int A = 0; A < s.nNodes(); A++ ) {
      for ( int i = 0; i < nDim; i++ )
        Y( i, A ) = s.grid.centers[s.assignedGridIndices[A]][i];
      dQ.template segment< nDim >( A * nBlock ) = G * Y.col( A );
    }
    Eigen::Matrix< double, nDim, 1 > Yc;
    s.particle->getCenterCoordinates( Yc.data() );

    s.particle->updatePhysicsExplicit( dQ.data(), dT, dT );
    Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() );
    s.particle->computeLumpedMomentum( P.data() );

    Eigen::Matrix< double, nDim, 1 >    momentum    = Eigen::Matrix< double, nDim, 1 >::Zero();
    Eigen::Matrix< double, nDim, nDim > firstMoment = Eigen::Matrix< double, nDim, nDim >::Zero();
    for ( int A = 0; A < s.nNodes(); A++ ) {
      momentum += P.template segment< nDim >( A * nBlock );
      firstMoment += P.template segment< nDim >( A * nBlock ) * Y.col( A ).transpose();
    }

    const Eigen::Matrix< double, nDim, 1 > vc                     = G * Yc / dT;
    const double                           edge                   = 0.6;
    const double secondMomentOfTheDomain                          = withNSNIStabilization ? V0 * edge * edge / 12 : 0.0;
    const Eigen::Matrix< double, nDim, nDim > expectedFirstMoment = density * V0 * vc * Yc.transpose() +
                                                                    density * G / dT * secondMomentOfTheDomain;
    const double scale = density * V0 * vc.norm() * Yc.norm();

    throwExceptionOnFailure( ( momentum - density * V0 * vc ).norm() <= 1e-10 * density * V0 * vc.norm(),
                             name + ": the lumped momentum of a linear increment is not the particle's momentum" );
    throwExceptionOnFailure( ( firstMoment - expectedFirstMoment ).norm() <= 1e-10 * scale,
                             Marmot::Utilities::MakeString()
                               << name << ": the first moment of the lumped momentum differs by "
                               << ( firstMoment - expectedFirstMoment ).norm() << " of " << scale );
  }

} // namespace

void testUnknownParticleIsRejected()
{
  Marmot::Discretizations::Meshfree::ReproducingKernelApproximation approximation( 2, 1 );
  const auto                                                        v = vertices( "Quad" );
  throwExceptionOnFailure( throws( [&]() {
                             std::unique_ptr< Marmot::Discretizations::Meshfree::Particle > p(
                               Marmot::Registry::ParticleFactory::createParticle( "GradientEnhancedFiniteStrain/"
                                                                                  "NoSuchShape",
                                                                                  1,
                                                                                  v.data(),
                                                                                  static_cast< int >( v.size() ),
                                                                                  0.36,
                                                                                  "GRADIENTENHANCEDCOMPRESSIBLE"
                                                                                  "NEOHOOKEDAMAGE",
                                                                                  matProps.data(),
                                                                                  static_cast< int >( matProps.size() ),
                                                                                  approximation ) );
                           } ),
                           "an unknown particle type must throw" );
}

int main()
{
  const std::vector< std::pair< std::string, std::string > > planeStrain = {
    { "GradientEnhancedFiniteStrain/PlaneStrain/Point", "Point" },
    { "GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSNNI/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSQCNI_R/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSQCNI_RU/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSNNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSQCNI_RxNSNI/PlaneStrain/Quad", "Quad" },
    { "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/PlaneStrain/Quad", "Quad" },
  };
  const std::vector< std::pair< std::string, std::string > > solid = {
    { "GradientEnhancedFiniteStrain/3D/Point", "Point3D" },
    { "GradientEnhancedFiniteStrainSQCNI/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSNNI/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSQCNI_R/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSQCNI_RU/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSNNIxNSNI/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSQCNI_RxNSNI/3D/Hexa", "Hexa" },
    { "GradientEnhancedFiniteStrainSQCNI_RUxNSNI/3D/Hexa", "Hexa" },
  };

  std::vector< std::function< void() > > testFunctions;
  for ( const auto& [name, shape] : planeStrain ) {
    testFunctions.push_back( [name = name, shape = shape]() { checkParticle< 2 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkTangentWithdL_dN< 2 >( name, shape ); } );
    if ( shape == "Quad" )
      testFunctions.push_back( [name = name, shape = shape]() { checkWeakFormCorrection< 2 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 2 >( name, shape ); } );
  }
  for ( const auto& [name, shape] : solid ) {
    testFunctions.push_back( [name = name, shape = shape]() { checkParticle< 3 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkTangentWithdL_dN< 3 >( name, shape ); } );
    if ( shape == "Hexa" )
      testFunctions.push_back( [name = name, shape = shape]() { checkWeakFormCorrection< 3 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 3 >( name, shape ); } );
  }
  for ( const auto& [name, shape] : planeStrain )
    testFunctions.push_back(
      [name = name, shape = shape]() { checkOwnPointsMoveWithTheLinearizedInterpolation< 2 >( name, shape ); } );
  for ( const auto& [name, shape] : solid )
    testFunctions.push_back(
      [name = name, shape = shape]() { checkOwnPointsMoveWithTheLinearizedInterpolation< 3 >( name, shape ); } );
  testFunctions.push_back( testSQCNIxNSNIFirstIncrementIsDeterministic );

  // the explicit scheme, with one nonlocal variable ...
  testFunctions.push_back( []() {
    checkExplicitTotalFormAfterAnAcceptedIncrement< 2 >( "GradientEnhancedFiniteStrain/PlaneStrain/Point", "Point" );
    checkExplicitTotalFormAfterAnAcceptedIncrement< 3 >( "GradientEnhancedFiniteStrain/3D/Point", "Point3D" );
  } );
  for ( const auto& [name, shape] : planeStrain ) {
    testFunctions.push_back( [name = name, shape = shape]() {
      checkExplicitAgreesWithImplicit< 2, 1 >( name, shape, testMaterialName, matProps, 0.02 );
      checkLumpedMomentumOfALinearIncrement< 2, 1 >( name,
                                                     shape,
                                                     viscousMaterialName,
                                                     matProps,
                                                     matProps[5],
                                                     name.find( "NSNI" ) != std::string::npos );
      checkLumpedOperators< 2, 1 >( name,
                                    shape,
                                    viscousMaterialName,
                                    matProps,
                                    matProps[5],
                                    ViscousNonlocalField::eta,
                                    ViscousNonlocalField::m );
    } );
  }
  for ( const auto& [name, shape] : solid ) {
    testFunctions.push_back( [name = name, shape = shape]() {
      checkExplicitAgreesWithImplicit< 3, 1 >( name, shape, testMaterialName, matProps, 0.02 );
      checkLumpedMomentumOfALinearIncrement< 3, 1 >( name,
                                                     shape,
                                                     viscousMaterialName,
                                                     matProps,
                                                     matProps[5],
                                                     name.find( "NSNI" ) != std::string::npos );
      checkLumpedOperators< 3, 1 >( name,
                                    shape,
                                    viscousMaterialName,
                                    matProps,
                                    matProps[5],
                                    ViscousNonlocalField::eta,
                                    ViscousNonlocalField::m );
    } );
  }
  // ... and with six, GM7 in the frame of the polar rotation
  for ( const auto& [name, shape] :
        std::vector< std::pair< std::string, std::string > >{ { "GradientEnhanced6FiniteStrain/3D/Point", "Point3D" },
                                                              { "GradientEnhanced6FiniteStrainSQCNI/3D/Hexa", "Hexa" },
                                                              { "GradientEnhanced6FiniteStrainSQCNIxNSNI/3D/Hexa",
                                                                "Hexa" } } ) {
    testFunctions.push_back( [name = name, shape = shape]() {
      checkExplicitAgreesWithImplicit< 3, 6 >( name, shape, "GM7/COROTATIONAL", gm7Props, 2e-5 );
      checkLumpedOperators< 3, 6 >( name,
                                    shape,
                                    "GM7/COROTATIONAL",
                                    gm7Props,
                                    gm7Props[14],
                                    gm7Props[15],
                                    gm7Props[16] );
    } );
  }

  // the bound of the absolute stiffness, for every particle with an own inertia
  for ( const auto& [name, shape] : planeStrain )
    testFunctions.push_back( [name = name, shape = shape]() {
      checkBoundOfAbsoluteStiffness< 2, 1 >( name, shape, testMaterialName, matProps, 0.02, 1.5 );
    } );
  for ( const auto& [name, shape] : solid )
    testFunctions.push_back( [name = name, shape = shape]() {
      checkBoundOfAbsoluteStiffness< 3, 1 >( name, shape, testMaterialName, matProps, 0.02, 1.5 );
    } );
  for ( const auto& [name, shape] :
        std::vector< std::pair< std::string, std::string > >{ { "GradientEnhanced6FiniteStrain/3D/Point", "Point3D" },
                                                              { "GradientEnhanced6FiniteStrainSQCNI/3D/Hexa", "Hexa" },
                                                              { "GradientEnhanced6FiniteStrainSQCNIxNSNI/3D/Hexa",
                                                                "Hexa" } } )
    testFunctions.push_back( [name = name, shape = shape]() {
      checkBoundOfAbsoluteStiffness< 3, 6 >( name, shape, "GM7/COROTATIONAL", gm7Props, 2e-4, 1.5 );
    } );

  testFunctions.push_back( testUnknownParticleIsRejected );

  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
