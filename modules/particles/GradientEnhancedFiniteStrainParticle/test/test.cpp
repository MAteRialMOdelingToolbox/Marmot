#include "Marmot/GradientEnhancedFiniteStrainParticle.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Marmot::Meshfree;
using namespace Marmot::Testing;

namespace {

  // GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE: K, G, kappa0, kappaF, l, rho
  const std::vector< double > matProps = { 3500., 1500., 1e-3, 1e-2, 0.3, 2.0 };

  constexpr int nDim = 2;

  // a 5 x 5 grid of 2nd order B-spline kernels with spacing 1, and a particle in its interior
  struct KernelGrid {
    std::vector< std::array< double, 2 > >                                             centers;
    std::vector< std::unique_ptr< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed > > kernels;
    std::vector< const MarmotMeshfreeKernelFunction* >                                 pointers;

    KernelGrid()
    {
      for ( int j = 0; j < 5; j++ )
        for ( int i = 0; i < 5; i++ )
          centers.push_back( { double( i ), double( j ) } );
      for ( auto& c : centers ) {
        kernels.emplace_back(
          std::make_unique< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed >( c.data(), nDim, 1.6 ) );
        pointers.push_back( kernels.back().get() );
      }
    }

    // the kernels whose support covers the given point
    std::vector< const MarmotMeshfreeKernelFunction* > covering( const double* x ) const
    {
      std::vector< const MarmotMeshfreeKernelFunction* > result;
      for ( const auto* k : pointers )
        if ( k->isInSupport( x ) )
          result.push_back( k );
      return result;
    }
  };

  struct Setup {
    KernelGrid                                   grid;
    MarmotMeshfreeReproducingKernelApproximation approximation{ nDim, 1 };
    std::unique_ptr< MarmotParticle >            particle;
    std::vector< double >                        stateVars;
    int                                          nNodes = 0;

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

    void assignKernels()
    {
      // the union over all evaluation points, as the host's particle manager does it
      const int             nEval = particle->getNumberOfEvaluationPoints();
      std::vector< double > eval( nEval * nDim );
      particle->getEvaluationCoordinates( eval.data() );

      std::vector< const MarmotMeshfreeKernelFunction* > assigned;
      for ( const auto* k : grid.pointers )
        for ( int e = 0; e < nEval; e++ )
          if ( k->isInSupport( &eval[e * nDim] ) ) {
            assigned.push_back( k );
            break;
          }

      particle->assignMeshfreeKernelFunctions( assigned );
      nNodes = assigned.size();
    }

    int nDof() const { return nNodes * ( nDim + 1 ); }

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

  // interleaved particle dof vector [u_1 N_1 u_2 N_2 ...]
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

  const std::vector< double > pointVertex  = { 1.9, 2.2 };
  const std::vector< double > quadVertices = { 1.6, 1.9, 2.2, 1.9, 2.2, 2.5, 1.6, 2.5 };

  // exact: the full tangent; approximate (NSNI): the stabilization's tangent omits the derivative of
  // dTau/dDeltaF itself (d2tau/dF2 and the damage dependence), so only the N rows are exact there
  enum class Tangent { Exact, NSNIApproximate };

  void checkTangent( const std::string& what, const Eigen::MatrixXd& K, const Eigen::MatrixXd& numK, Tangent mode )
  {
    const auto be = blockErrors( K, numK );
    if ( mode == Tangent::Exact ) {
      const double err = ( K - numK ).norm() / numK.norm();
      throwExceptionOnFailure( err < 1e-7, MakeString() << what << ": tangent inconsistent, relative error " << err );
    }
    else {
      throwExceptionOnFailure( be[2] < 1e-5 && be[3] < 1e-7,
                               MakeString() << what << ": nonlocal rows inconsistent, " << be[2] << ", " << be[3] );
      throwExceptionOnFailure( be[0] < 5e-3,
                               MakeString() << what << ": U-U block beyond the NSNI approximation, " << be[0] );
    }
  }

  void checkParticle( const std::string& name, const std::vector< double >& vertices, double volume, Tangent mode )
  {
    // first step, from the undeformed state
    {
      Setup                 s( name, vertices, volume );
      const Eigen::VectorXd dQ = increment( s.nNodes, 0.02, 3e-3 );
      const auto [P, K]        = s.trial( dQ );
      const auto numK          = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
      checkTangent( name + ", first step", K, numK, mode );

      // translation invariance: the displacement residuals of every direction sum to zero
      for ( int i = 0; i < nDim; i++ ) {
        double sum = 0;
        for ( int A = 0; A < s.nNodes; A++ )
          sum += P[A * ( nDim + 1 ) + i];
        throwExceptionOnFailure( std::abs( sum ) < 1e-10 * P.norm(), name + ": residuals are not self-equilibrated" );
      }
    }

    // second step, from an accepted finite deformation with damage history
    {
      Setup s( name, vertices, volume );
      s.accept( increment( s.nNodes, 0.02, 3e-3 ) );
      const Eigen::VectorXd dQ = increment( s.nNodes, 0.01, 1e-3 );
      const auto [P, K]        = s.trial( dQ );
      const auto numK          = numericalTangent( [&]( const Eigen::VectorXd& q ) { return s.trial( q ).first; }, dQ );
      checkTangent( name + ", second step", K, numK, mode );
    }

    // a finite rigid rotation of the kernel nodes does not load the particle
    {
      Setup           s( name, vertices, volume );
      const double    phi = 0.5;
      Eigen::Matrix2d R;
      R << std::cos( phi ), -std::sin( phi ), std::sin( phi ), std::cos( phi );

      // nodal coefficients reproducing u = (R - I) X exactly (linear completeness)
      Eigen::VectorXd dQ = Eigen::VectorXd::Zero( s.nDof() );
      int             A  = 0;
      // the assigned kernels are in grid order; recover their centres through the grid
      std::vector< double > eval( s.particle->getNumberOfEvaluationPoints() * nDim );
      s.particle->getEvaluationCoordinates( eval.data() );
      for ( size_t g = 0; g < s.grid.pointers.size(); g++ ) {
        bool covers = false;
        for ( size_t e = 0; e < eval.size() / nDim; e++ )
          covers = covers || s.grid.pointers[g]->isInSupport( &eval[e * nDim] );
        if ( !covers )
          continue;
        const Eigen::Vector2d X( s.grid.centers[g].data() );
        dQ.segment< 2 >( A * ( nDim + 1 ) ) = ( R - Eigen::Matrix2d::Identity() ) * X;
        A++;
      }

      const auto P = s.trial( dQ ).first;
      throwExceptionOnFailure( P.norm() < 1e-8,
                               MakeString() << name << ": rigid rotation loads the particle, |P| = " << P.norm() );
    }
  }

} // namespace

void testPointParticle()
{
  checkParticle( "GradientEnhancedFiniteStrain/PlaneStrain/Point", pointVertex, 0.36, Tangent::Exact );
}

void testSQCNIParticle()
{
  checkParticle( "GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad", quadVertices, 0.36, Tangent::Exact );
}

void testSQCNIxNSNIParticle()
{
  checkParticle( "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad",
                 quadVertices,
                 0.36,
                 Tangent::NSNIApproximate );
}

void testSQCNIxNSNIFirstIncrementIsDeterministic()
{
  // regression: the intermediate second moments used to be set only on the first accepted increment, so the
  // stabilization of the very first increment was scaled by whatever the heap held. Two particles built after
  // different heap histories must agree bit for bit.
  const std::string name = "GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad";

  Setup      a( name, quadVertices, 0.36 );
  const auto Pa = a.trial( increment( a.nNodes, 0.02, 3e-3 ) ).first;

  {
    std::vector< double > garbage( 1 << 16, 1.2345e300 ); // dirty the heap, then free it for reuse
  }
  Setup      b( name, quadVertices, 0.36 );
  const auto Pb = b.trial( increment( b.nNodes, 0.02, 3e-3 ) ).first;

  throwExceptionOnFailure( Pa.allFinite() && ( Pa - Pb ).norm() == 0.0,
                           "first-increment residual depends on uninitialized memory" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testPointParticle,
                                                               testSQCNIParticle,
                                                               testSQCNIxNSNIParticle,
                                                               testSQCNIxNSNIFirstIncrementIsDeterministic };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
