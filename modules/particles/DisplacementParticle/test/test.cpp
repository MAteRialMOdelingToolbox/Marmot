#include "Marmot/DisplacementParticle.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
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
    int nDof() const { return nNodes() * nDim; }

    std::pair< Eigen::VectorXd, Eigen::MatrixXd > trial( const Eigen::VectorXd& dQ )
    {
      const auto      backup = stateVars;
      Eigen::VectorXd P      = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K      = Eigen::MatrixXd::Zero( nDof(), nDof() );
      particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
      stateVars = backup;
      return { P, K };
    }

    void accept( const Eigen::VectorXd& dQ )
    {
      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof() );
      Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof(), nDof() );
      particle->computePhysicsKernels( dQ.data(), P.data(), K.data(), 1.0, 1.0 );
      particle->acceptStateAndPosition();
      assignKernels();
    }
  };

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

  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
