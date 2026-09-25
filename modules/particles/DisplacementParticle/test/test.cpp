#include "Marmot/DisplacementParticle.h"
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

  // COMPRESSIBLENEOHOOKE: K, G, rho
  const std::vector< double > matProps = { 3500., 1500., 2.0 };

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
    KernelGrid< nDim >                           grid;
    MarmotMeshfreeReproducingKernelApproximation approximation{ nDim, 1 };
    std::unique_ptr< MarmotParticle >            particle;
    std::vector< double >                        stateVars;
    std::vector< int >                           assignedGridIndices;
    double                                       volume;

    Setup( const std::string& name, const std::vector< double >& vertices, double volume_ ) : volume( volume_ )
    {
      particle.reset( MarmotLibrary::MarmotParticleFactory::createParticle( name,
                                                                            1,
                                                                            vertices.data(),
                                                                            vertices.size(),
                                                                            volume,
                                                                            "COMPRESSIBLENEOHOOKE",
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

    // a uniform translation t of all nodes
    Eigen::VectorXd translation( const Eigen::Matrix< double, nDim, 1 >& t ) const
    {
      Eigen::VectorXd dQ( nDof() );
      for ( int A = 0; A < nNodes(); A++ )
        dQ.template segment< nDim >( A * nDim ) = t;
      return dQ;
    }
    int nDof() const { return nNodes() * nDim; }

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
  };

  // fill freed heap chunks of many sizes with garbage, so that members which are read before they are set show up as
  // garbage instead of as zeros or as stale values of a previous, identical object
  void dirtyHeap()
  {
    for ( int n = 1; n <= 1024; n += 1 + n / 8 )
      std::vector< double >( n, 1.2345e300 );
  }

  Eigen::VectorXd increment( int nDof, double scale )
  {
    Eigen::VectorXd dQ( nDof );
    for ( int i = 0; i < nDof; i++ )
      dQ[i] = scale * std::sin( 1.7 * ( i + 1 ) );
    return dQ;
  }

  Eigen::MatrixXd numericalTangent( const std::function< Eigen::VectorXd( const Eigen::VectorXd& ) >& f,
                                    const Eigen::VectorXd&                                            Q0 )
  {
    const double    h = 1e-7;
    Eigen::MatrixXd numK( f( Q0 ).size(), Q0.size() );
    for ( int j = 0; j < Q0.size(); j++ ) {
      Eigen::VectorXd Qp = Q0, Qm = Q0;
      Qp[j] += h;
      Qm[j] -= h;
      numK.col( j ) = ( f( Qp ) - f( Qm ) ) / ( 2 * h );
    }
    return numK;
  }

  // the particle geometries: a point, a quad and a hexahedron in the interior of the kernel grid
  std::vector< double > vertices( const std::string& shape )
  {
    if ( shape == "Point" )
      return { 1.9, 2.2 };
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

  // faceted particles compute their volume from the vertices and require 0 here
  double volumeArgument( const std::string& shape )
  {
    return shape == "Point" ? 0.36 : 0.0;
  }

  template < int nDim >
  void checkParticle( const std::string& name, const std::string& shape, double tangentTolerance )
  {
    const auto   v   = vertices( shape );
    const double vol = volumeArgument( shape );

    // consistent tangent, from the undeformed state and after an accepted finite increment
    for ( bool afterStep : { false, true } ) {
      Setup< nDim > s( name, v, vol );
      if ( afterStep )
        s.accept( increment( s.nDof(), 0.03 ) );
      const Eigen::VectorXd dQ = increment( s.nDof(), afterStep ? 0.01 : 0.03 );
      const auto [P, K]        = s.trial( dQ );
      const auto   numK        = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
      const double err         = ( K - numK ).norm() / numK.norm();
      throwExceptionOnFailure( err < tangentTolerance,
                               MakeString() << name << ( afterStep ? ", second" : ", first" )
                                            << " step: tangent inconsistent, relative error " << err );

      // translation invariance: the residuals of every direction sum to zero
      for ( int i = 0; i < nDim; i++ ) {
        double sum = 0;
        for ( int A = 0; A < s.nNodes(); A++ )
          sum += P[A * nDim + i];
        throwExceptionOnFailure( std::abs( sum ) < 1e-10 * P.norm(), name + ": residuals are not self-equilibrated" );
      }
    }

    // a finite rigid rotation (about z) does not load the particle
    {
      Setup< nDim >   s( name, v, vol );
      const double    phi = 0.5;
      Eigen::Matrix3d R   = Eigen::Matrix3d::Identity();
      R.topLeftCorner< 2, 2 >() << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );
      Eigen::VectorXd dQ( s.nDof() );
      for ( int A = 0; A < s.nNodes(); A++ ) {
        const Eigen::Vector3d X( s.grid.centers[s.assignedGridIndices[A]].data() );
        dQ.segment< nDim >( A * nDim ) = ( ( R - Eigen::Matrix3d::Identity() ) * X ).template head< nDim >();
      }
      const auto P = s.trial( dQ ).first;
      throwExceptionOnFailure( P.norm() < 1e-8,
                               MakeString() << name << ": rigid rotation loads the particle, |P| = " << P.norm() );
    }

    // lumped mass (lumped inertia is optional in the particle interface). Body loads are not checked: the
    // particles advertise BODYFORCE, but computeBodyLoad is not implemented (empty, or throws for the SDI
    // variants), and EdelweissMeshfree applies body loads to cells only.
    {
      dirtyHeap(); // the density must be known right after initialization, before the first increment
      Setup< nDim >   s( name, v, vol );
      const double    V           = s.particle->getVolumeUndeformed();
      Eigen::VectorXd m           = Eigen::VectorXd::Zero( s.nDof() );
      bool            implemented = true;
      try {
        s.particle->computeLumpedInertia( m.data() );
      }
      catch ( const std::runtime_error& ) {
        implemented = false;
      }
      if ( implemented )
        throwExceptionOnFailure( std::abs( m.sum() - nDim * matProps[2] * V ) < 1e-12 * ( 1. + V ),
                                 MakeString() << name << ": lumped mass " << m.sum() << " != nDim rho V" );
    }

    // pressure load tangent on the first face of particles that have faces
    if ( shape != "Point" ) {
      Setup< nDim >         s( name, v, vol );
      const auto            it   = s.particle->getSupportedDistributedLoadTypes().find( "PRESSURE" );
      const double          p    = 2.0;
      const Eigen::VectorXd dQ   = increment( s.nDof(), 0.03 );
      auto                  load = [&]( const Eigen::VectorXd& q ) {
        s.trial( q ); // no state is committed; loads read the current configuration from the kernels
        Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() );
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
        // the configuration of the current increment is the one the particle computed last
        s.particle->computePhysicsKernels( q.data(), P.data(), K.data(), 1.0, 1.0 );
        P.setZero();
        K.setZero();
        s.particle->computeDistributedLoad( it->second, 1, &p, P.data(), K.data(), 1.0, 1.0 );
        return std::make_pair( P, K );
      };
      const auto backup   = s.stateVars;
      const auto [P0, K0] = load( dQ );
      const auto numK     = numericalTangent(
        [&]( const Eigen::VectorXd& q ) {
          const auto r = load( q ).first;
          s.stateVars  = backup;
          return r;
        },
        dQ );
      s.stateVars = backup;
      throwExceptionOnFailure( P0.norm() > 0, name + ": pressure must load the particle" );
      const double err = std::min( ( K0 - numK ).norm(), ( K0 + numK ).norm() ) / numK.norm();
      throwExceptionOnFailure( err < 1e-6, MakeString() << name << ": pressure load tangent inconsistent, " << err );
    }
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

  int numberOfFaces( const std::string& shape )
  {
    return shape == "Quad" ? 4 : shape == "Hexa" ? 6 : 0;
  }

  // the parts of the particle interface beyond the mechanics: properties, geometry, state, VCI, explicit kernels
  template < int nDim >
  void checkInterface( const std::string& name, const std::string& shape )
  {
    using Vec      = Eigen::Matrix< double, nDim, 1 >;
    const auto   v = vertices( shape );
    const double V = volumeArgument( shape );

    // properties
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

      // the number of VCI constraints is the size of the complete basis in nDim dimensions
      s.setProperties( { { "VCI order", 1.0 } } );
      throwExceptionOnFailure( s.particle->vci_getNumberOfConstraints() == nDim + 1,
                               MakeString() << name << ": " << s.particle->vci_getNumberOfConstraints()
                                            << " linear VCI constraints in " << nDim << "D" );

      throwExceptionOnFailure( s.particle->getFields() == std::vector< std::string >{ "displacement" },
                               name + ": fields" );
      throwExceptionOnFailure( s.particle->getNBaseDof() == nDim && s.particle->getDimension() == nDim,
                               name + ": dimension" );
      throwExceptionOnFailure( !s.particle->getParticleShape().empty(), name + ": shape" );
      throwExceptionOnFailure( throws( [&]() { s.particle->setInitialCondition( "no such condition", &one ); } ),
                               name + ": an unknown initial condition must throw" );
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

    // a translation is carried by every vertex and by the material point(s)
    {
      Setup< nDim > s( name, v, V );
      const Vec     t = Vec::LinSpaced( 0.1, 0.1 * nDim );
      s.accept( s.translation( t ) );

      const auto u = s.particle->getStateView( "displacement", 0 );
      throwExceptionOnFailure( ( Eigen::Map< const Vec >( u.stateLocation ) - t ).norm() < 1e-12,
                               name + ": displacement of a translation" );
      if ( numberOfFaces( shape ) > 0 )
        for ( const std::string& state : { "vertex displacements", "smoothing vertex displacements" } ) {
          const auto view = s.particle->getStateView( state, 0 );
          throwExceptionOnFailure( view.stateSize == int( v.size() ), name + ": size of " + state );
          for ( size_t k = 0; k < v.size() / nDim; k++ )
            throwExceptionOnFailure( ( Eigen::Map< const Vec >( view.stateLocation + k * nDim ) - t ).norm() < 1e-12,
                                     name + ": " + state + " of a translation" );
        }

      Vec translated;
      s.particle->getCenterCoordinates( translated.data() );
      Vec mean = Vec::Zero();
      for ( size_t k = 0; k < v.size() / nDim; k++ )
        mean += Eigen::Map< const Vec >( &v[k * nDim] ) / ( v.size() / nDim );
      throwExceptionOnFailure( ( translated - mean - t ).norm() < 1e-12, name + ": translated center" );
    }

    // VCI of order 0: the correction eta = M^-1 R makes the corrected test functions satisfy the divergence
    // theorem, R = int_dOmega T n - int_Omega grad T = 0, for every node whose correction is active (M > 0)
    if ( numberOfFaces( shape ) > 0 ) {
      // the host corrects before the first increment: nothing may depend on a previously accepted state. Dirty the
      // heap first, so that a basis which is never evaluated shows up as garbage instead of as zeros.
      dirtyHeap();
      Setup< nDim > s( name, v, V );
      const int     nC = s.particle->vci_getNumberOfConstraints();
      throwExceptionOnFailure( nC == 1, name + ": one constraint of order 0" );

      auto residual = [&]() {
        Eigen::VectorXd R = Eigen::VectorXd::Zero( s.nDof() ), vol = Eigen::VectorXd::Zero( s.nDof() );
        for ( int f = 1; f <= numberOfFaces( shape ); f++ )
          s.particle->vci_compute_Test_P_BoundaryIntegral( R.data(), nullptr, f );
        s.particle->vci_compute_TestGradient_P_Integral( vol.data() );
        s.particle->vci_compute_Test_PGradient_Integral( vol.data() );
        return Eigen::VectorXd( R - vol );
      };

      Eigen::VectorXd M = Eigen::VectorXd::Zero( s.nNodes() );
      s.particle->vci_compute_MMatrix( M.data() );
      const Eigen::VectorXd R0 = residual();
      Eigen::VectorXd       eta( s.nDof() );
      for ( int A = 0; A < s.nNodes(); A++ )
        eta.template segment< nDim >( A * nDim ) = M[A] > 0
                                                     ? Eigen::VectorXd( R0.template segment< nDim >( A * nDim ) / M[A] )
                                                     : Eigen::VectorXd::Zero( nDim );
      throwExceptionOnFailure( std::abs( M.maxCoeff() - s.particle->getVolumeUndeformed() ) < 1e-12,
                               name + ": M of order 0 is the volume" );
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
      Eigen::VectorXd R = Eigen::VectorXd::Zero( s.nDof() );
      const Vec       n = Vec::LinSpaced( 0.3, 0.6 );
      s.particle->acceptStateAndPosition();
      s.particle->vci_compute_Test_P_BoundaryIntegral( R.data(), n.data(), 0 );
      Vec sum = Vec::Zero();
      for ( int A = 0; A < s.nNodes(); A++ )
        sum += R.template segment< nDim >( A * nDim );
      throwExceptionOnFailure( ( sum - n ).norm() < 1e-12, name + ": point boundary integral" );
    }

    // dynamics: with Newmark integration, the inertia enters residual and tangent consistently
    {
      Setup< nDim > s( name, v, V );
      s.setProperties( { { "newmark-beta beta", 0.25 }, { "newmark-beta gamma", 0.5 } } );
      const double          dT = 0.01;
      const Eigen::VectorXd dQ = increment( s.nDof(), 0.003 );
      const auto [P, K]        = s.trial( dQ, dT );
      const auto   numK = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q, dT ).first; }, dQ );
      const double err  = ( K - numK ).norm() / numK.norm();
      throwExceptionOnFailure( err < ( name.find( "NSNI" ) != std::string::npos ? 5e-3 : 1e-7 ),
                               MakeString() << name << ": dynamic tangent inconsistent, " << err );

      // lumped momentum of the accepted state: m v
      s.accept( dQ, dT );
      Eigen::VectorXd p = Eigen::VectorXd::Zero( s.nDof() );
      if ( !throws( [&]() { s.particle->computeLumpedMomentum( p.data() ); } ) ) {
        const auto vel = s.particle->getStateView( "velocity", 0 );
        for ( int i = 0; i < nDim; i++ ) {
          double sum = 0;
          for ( int A = 0; A < s.nNodes(); A++ )
            sum += p[A * nDim + i];
          throwExceptionOnFailure( std::abs( sum - matProps[2] * s.particle->getVolumeUndeformed() *
                                                     vel.stateLocation[i] ) < 1e-10 * ( 1 + std::abs( sum ) ),
                                   name + ": lumped momentum is not m v" );
        }
      }
    }

    // the explicit internal force (NSNI particles implement the explicit interface) is the quasi-static implicit one
    if ( name.find( "NSNI" ) != std::string::npos ) {
      Setup< nDim >         s( name, v, V );
      const Eigen::VectorXd dQ        = increment( s.nDof(), 0.03 );
      const Eigen::VectorXd P         = s.trial( dQ ).first;
      Eigen::VectorXd       PExplicit = Eigen::VectorXd::Zero( s.nDof() );
      s.particle->updatePhysicsExplicit( dQ.data(), 1.0, 1.0 );
      s.particle->computePhysicsKernelsExplicit( PExplicit.data() );
      throwExceptionOnFailure( ( P - PExplicit ).norm() < 1e-10 * P.norm(),
                               MakeString() << name << ": explicit internal force differs from the implicit one, "
                                            << ( P - PExplicit ).norm() / P.norm() );
    }

    // the weak-form correction on all faces balances the internal force of a homogeneous deformation (divergence
    // theorem of the smoothed gradient), before and after an accepted finite increment
    // (after a step, only the full SQCNI smooths over the deformed particle itself, which makes the balance exact)
    const bool smoothsOverDeformedParticle = name.find( "SQCNI/" ) != std::string::npos ||
                                             name.find( "SQCNIxNSNI/" ) != std::string::npos;
    if ( numberOfFaces( shape ) > 0 &&
         Setup< nDim >( name, v, V ).particle->getSupportedDistributedLoadTypes().count( "CWFCORRECTION" ) )
      for ( bool afterStep : { false, true } ) {
        if ( afterStep && !smoothsOverDeformedParticle )
          continue;
        Setup< nDim >                       s( name, v, V );
        Eigen::Matrix< double, nDim, nDim > G = Eigen::Matrix< double, nDim, nDim >::Zero();
        G.diagonal().setConstant( 0.04 );
        G( 0, 1 )        = 0.03;
        auto homogeneous = [&]( const Eigen::Matrix< double, nDim, nDim >& H ) {
          Eigen::VectorXd dQ( s.nDof() );
          for ( int A = 0; A < s.nNodes(); A++ )
            dQ.template segment< nDim >( A * nDim ) = H * Eigen::Map< const Vec >(
                                                            s.grid.centers[s.assignedGridIndices[A]].data() );
          return dQ;
        };
        if ( afterStep )
          s.accept( homogeneous( G ) );
        const Eigen::VectorXd dQ = homogeneous( 0.5 * G.transpose() );
        // the loads read the configuration of the last computation, so compute once more without restoring
        Eigen::VectorXd P = Eigen::VectorXd::Zero( s.nDof() ), Pc = P;
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() ), Kc = K;
        s.particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
        const int cwf = s.particle->getSupportedDistributedLoadTypes().at( "CWFCORRECTION" );
        for ( int f = 1; f <= numberOfFaces( shape ); f++ )
          s.particle->computeDistributedLoad( cwf, f, nullptr, Pc.data(), Kc.data(), 1.0, 1.0 );
        throwExceptionOnFailure( ( P + Pc ).norm() < 1e-10 * P.norm(),
                                 MakeString() << name << ( afterStep ? ", second step" : ", first step" )
                                              << ": weak-form correction does not balance a homogeneous state, "
                                              << ( P + Pc ).norm() / P.norm() );
      }

    // distributed loads in the explicit interface agree with the implicit ones
    if ( numberOfFaces( shape ) > 0 ) {
      const auto& loads = Setup< nDim >( name, v, V ).particle->getSupportedDistributedLoadTypes();
      for ( const auto& [loadName, type] : loads ) {
        Setup< nDim >         s( name, v, V );
        const double          p  = 2.0;
        const Eigen::VectorXd dQ = increment( s.nDof(), 0.03 );
        Eigen::VectorXd       P = Eigen::VectorXd::Zero( s.nDof() ), PExplicit = P;
        Eigen::MatrixXd       K = Eigen::MatrixXd::Zero( s.nDof(), s.nDof() );
        s.trial( dQ );
        s.particle->computeDistributedLoad( type, 2, &p, P.data(), K.data(), 1.0, 1.0 );
        // only the pressure is available explicitly, and the SDI particles inherit the empty explicit interface
        if ( throws(
               [&]() { s.particle->computeDistributedLoadExplicit( type, 2, &p, PExplicit.data(), 1.0, 1.0 ); } ) ||
             PExplicit.norm() == 0 )
          continue;
        throwExceptionOnFailure( ( P - PExplicit ).norm() < 1e-12 * ( 1 + P.norm() ),
                                 name + ": explicit " + loadName + " differs from the implicit one" );
      }
    }
  }

} // namespace

int main()
{
  const std::vector< std::pair< std::string, std::string > > planeStrain = {
    { "Displacement/PlaneStrain/Point", "Point" },
    { "DisplacementSQCNI/PlaneStrain/Quad", "Quad" },
    { "DisplacementSNNI/PlaneStrain/Quad", "Quad" },
    { "DisplacementSQCNI_R/PlaneStrain/Quad", "Quad" },
    { "DisplacementSQCNI_RU/PlaneStrain/Quad", "Quad" },
    { "Displacement/SQCNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "Displacement/SNNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "Displacement/R-SNNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "Displacement/RS-SNNIxNSNI/PlaneStrain/Quad", "Quad" },
    { "DisplacementSQCNIxSDI/PlaneStrain/Quad", "Quad" },
  };
  const std::vector< std::pair< std::string, std::string > > solid = {
    { "DisplacementSQCNI/3D/Hexa", "Hexa" },
    { "Displacement/SQCNIxNSNI/3D/Hexa", "Hexa" },
    { "Displacement/SNNIxNSNI/3D/Hexa", "Hexa" },
    { "Displacement/R-SNNIxNSNI/3D/Hexa", "Hexa" },
    { "Displacement/RS-SNNIxNSNI/3D/Hexa", "Hexa" },
    { "DisplacementSQCNIxSDI/3D/Hexa", "Hexa" },
    { "DisplacementSNNIxSDI/3D/Hexa", "Hexa" },
    { "DisplacementR-SNNIxSDI/3D/Hexa", "Hexa" },
    { "DisplacementRS-SNNIxSDI/3D/Hexa", "Hexa" },
  };

  std::vector< std::function< void() > > testFunctions;
  // the NSNI stabilization's tangent omits d2tau/dF2 (not exposed by the material interface), so it is
  // approximate by construction: measured 1.1e-3 (2D) and 1.3e-3 (3D) relative error; all others are exact
  auto tolerance = []( const std::string& name ) { return name.find( "NSNI" ) != std::string::npos ? 5e-3 : 1e-7; };
  for ( const auto& [name, shape] : planeStrain )
    testFunctions.push_back(
      [&, name = name, shape = shape]() { checkParticle< 2 >( name, shape, tolerance( name ) ); } );
  for ( const auto& [name, shape] : solid )
    testFunctions.push_back(
      [&, name = name, shape = shape]() { checkParticle< 3 >( name, shape, tolerance( name ) ); } );
  for ( const auto& [name, shape] : planeStrain )
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 2 >( name, shape ); } );
  for ( const auto& [name, shape] : solid )
    testFunctions.push_back( [name = name, shape = shape]() { checkInterface< 3 >( name, shape ); } );

  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
