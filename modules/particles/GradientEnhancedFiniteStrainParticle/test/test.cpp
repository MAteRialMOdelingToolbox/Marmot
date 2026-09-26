#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Meshfree;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };

  // a grid of 5^nDim 2nd order B-spline kernels with spacing 1 around the particle
  template < int nDim >
  struct KernelGrid {
    std::vector< std::array< double, 3 > >                                             centers;
    std::vector< std::unique_ptr< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed > > kernels;
    std::vector< const MarmotMeshfreeKernelFunction* >                                 pointers;

    KernelGrid()
    {
      const int nz = nDim == 3 ? 5 : 1;
      for ( int k = 0; k < nz; k++ )
        for ( int j = 0; j < 5; j++ )
          for ( int i = 0; i < 5; i++ )
            centers.push_back( { double( i ), double( j ), double( k ) } );
      for ( auto& c : centers ) {
        kernels.emplace_back(
          std::make_unique< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed >( c.data(), nDim, 1.6 ) );
        pointers.push_back( kernels.back().get() );
      }
    }
  };

  template < int nDim >
  struct Setup {
    static constexpr int nBlock = nDim + 1; // interleaved particle dof vector [u_1 N_1 u_2 N_2 ...]

    KernelGrid< nDim >                           grid;
    MarmotMeshfreeReproducingKernelApproximation approximation{ nDim, 1 };
    std::unique_ptr< MarmotParticle >            particle;
    std::vector< double >                        stateVars;
    std::vector< int >                           assignedGridIndices;

    Setup( const std::string& name, const std::vector< double >& vertices, double volume )
    {
      particle.reset(
        MarmotLibrary::MarmotParticleFactory::createParticle( name,
                                                              1,
                                                              vertices.data(),
                                                              vertices.size(),
                                                              volume,
                                                              "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                                              matProps.data(),
                                                              matProps.size(),
                                                              approximation ) );
      stateVars.assign( particle->getNumberOfRequiredStateVars(), 0.0 );
      particle->assignStateVars( stateVars.data(), stateVars.size() );
      particle->initializeYourself();
      assignKernels();
    }

    // the union over all evaluation points, as the host's particle manager does it
    void assignKernels()
    {
      const int             nEval = particle->getNumberOfEvaluationPoints();
      std::vector< double > eval( nEval * nDim );
      particle->getEvaluationCoordinates( eval.data() );

      std::vector< const MarmotMeshfreeKernelFunction* > assigned;
      assignedGridIndices.clear();
      for ( size_t g = 0; g < grid.pointers.size(); g++ )
        for ( int e = 0; e < nEval; e++ )
          if ( grid.pointers[g]->isInSupport( &eval[e * nDim] ) ) {
            assigned.push_back( grid.pointers[g] );
            assignedGridIndices.push_back( g );
            break;
          }
      particle->assignMeshfreeKernelFunctions( assigned );
    }

    int nNodes() const { return assignedGridIndices.size(); }
    int nDof() const { return nNodes() * nBlock; }

    // set properties by name, the others to zero
    void setProperties( const std::map< std::string, double >& values )
    {
      const auto            names = particle->getPropertyNames();
      std::vector< double > p( names.size(), 0.0 );
      for ( size_t i = 0; i < names.size(); i++ )
        if ( values.count( names[i] ) )
          p[i] = values.at( names[i] );
      particle->setProperties( p.data(), p.size() );
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

  template < int nDim >
  Eigen::VectorXd increment( int nNodes, double scaleU, double N0 )
  {
    Eigen::VectorXd dQ( nNodes * ( nDim + 1 ) );
    for ( int A = 0; A < nNodes; A++ ) {
      for ( int i = 0; i < nDim; i++ )
        dQ[A * ( nDim + 1 ) + i] = scaleU * std::sin( 1.7 * ( A * nDim + i + 1 ) );
      dQ[A * ( nDim + 1 ) + nDim] = N0 * ( 1. + 0.3 * std::cos( 0.9 * ( A + 1 ) ) );
    }
    return dQ;
  }

  template < int nDim >
  Eigen::MatrixXd numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                    const Eigen::VectorXd&                                            Q0 )
  {
    Eigen::MatrixXd numK( f( Q0 ).size(), Q0.size() );
    for ( int j = 0; j < Q0.size(); j++ ) {
      const bool      isNonlocal = ( j % ( nDim + 1 ) ) == nDim;
      const double    h          = isNonlocal ? 1e-9 : 1e-7;
      Eigen::VectorXd Qp = Q0, Qm = Q0;
      Qp[j] += h;
      Qm[j] -= h;
      numK.col( j ) = ( f( Qp ) - f( Qm ) ) / ( 2 * h );
    }
    return numK;
  }

  // relative errors of the U-U, U-N, N-U, N-N blocks in the interleaved layout
  template < int nDim >
  std::array< double, 4 > blockErrors( const Eigen::MatrixXd& K, const Eigen::MatrixXd& numK )
  {
    std::array< double, 4 > num{}, den{};
    for ( int a = 0; a < K.rows(); a++ )
      for ( int b = 0; b < K.cols(); b++ ) {
        const int blk = 2 * ( a % ( nDim + 1 ) == nDim ) + ( b % ( nDim + 1 ) == nDim );
        num[blk] += std::pow( K( a, b ) - numK( a, b ), 2 );
        den[blk] += std::pow( numK( a, b ), 2 );
      }
    std::array< double, 4 > err;
    for ( int i = 0; i < 4; i++ )
      err[i] = std::sqrt( num[i] / std::max( den[i], 1e-300 ) );
    return err;
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
  void checkTangent( const std::string& what, const Eigen::MatrixXd& K, const Eigen::MatrixXd& numK, Tangent mode )
  {
    const auto be = blockErrors< nDim >( K, numK );
    if ( mode == Tangent::Exact ) {
      const double err = ( K - numK ).norm() / numK.norm();
      throwExceptionOnFailure( err < 1e-7, MakeString() << what << ": tangent inconsistent, relative error " << err );
    }
    else {
      throwExceptionOnFailure( be[2] < 5e-5 && be[3] < 1e-7,
                               MakeString() << what << ": nonlocal rows inconsistent, " << be[2] << ", " << be[3] );
      throwExceptionOnFailure( be[0] < 5e-3,
                               MakeString() << what << ": U-U block beyond the NSNI approximation, " << be[0] );
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
                               MakeString() << name << ": rigid rotation loads the particle, |P| = " << P.norm() );
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
        const double err = std::min( ( K0 - numK ).norm(), ( K0 + numK ).norm() ) / numK.norm();
        throwExceptionOnFailure( err < 1e-6,
                                 MakeString() << name << ": " << loadName << " load tangent inconsistent, " << err );
        for ( int A = 0; A < s.nNodes(); A++ )
          throwExceptionOnFailure( P0[A * nB + nDim] == 0,
                                   name + ": " + loadName + " must not load the nonlocal field" );
      }
    }
  }

  // the weak-form correction on all faces balances the displacement residual of a homogeneous deformation
  // (divergence theorem of the smoothed gradient); after a step, only the full SQCNI smooths over the deformed
  // particle itself, which makes the balance exact
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
                               MakeString() << name << ( afterStep ? ", second step" : ", first step" )
                                            << ": weak-form correction does not balance a homogeneous state, "
                                            << std::sqrt( res / ref ) );
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
      throwExceptionOnFailure( throws( [&]() { s.particle->setProperties( tooFew.data(), tooFew.size() ); } ),
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
                               MakeString()
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
                                   MakeString() << name << ": VCI residual after the correction, "
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

void testUnknownParticleIsRejected()
{
  MarmotMeshfreeReproducingKernelApproximation approximation( 2, 1 );
  const auto                                   v = vertices( "Quad" );
  throwExceptionOnFailure( throws( [&]() {
                             std::unique_ptr< MarmotParticle > p(
                               MarmotLibrary::MarmotParticleFactory::createParticle( "GradientEnhancedFiniteStrain/"
                                                                                     "NoSuchShape",
                                                                                     1,
                                                                                     v.data(),
                                                                                     v.size(),
                                                                                     0.36,
                                                                                     "GRADIENTENHANCEDCOMPRESSIBLE"
                                                                                     "NEOHOOKEDAMAGE",
                                                                                     matProps.data(),
                                                                                     matProps.size(),
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
    if ( shape == "Quad" )
      testFunctions.push_back( [name = name, shape = shape]() { checkWeakFormCorrection< 2 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 2 >( name, shape ); } );
  }
  for ( const auto& [name, shape] : solid ) {
    testFunctions.push_back( [name = name, shape = shape]() { checkParticle< 3 >( name, shape ); } );
    if ( shape == "Hexa" )
      testFunctions.push_back( [name = name, shape = shape]() { checkWeakFormCorrection< 3 >( name, shape ); } );
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 3 >( name, shape ); } );
  }
  testFunctions.push_back( testSQCNIxNSNIFirstIncrementIsDeterministic );
  testFunctions.push_back( testUnknownParticleIsRejected );

  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
