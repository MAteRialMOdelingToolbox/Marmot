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
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLagrangianCellGeometry.h"
#include "Marmot/MarmotMPMLibrary.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximationImplicit.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotParticleDomain.h"
#include "Marmot/MarmotParticleLibrary.h"
#include "Marmot/MarmotParticleStiffnessVisitors.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::Discretizations::Meshfree;
using namespace Marmot::Testing;

namespace {

  // all exponent tuples (a_1, ..., a_dim) with sum <= order
  std::vector< std::vector< int > > exponents( int order, int dim )
  {
    if ( dim == 0 )
      return { {} };
    std::vector< std::vector< int > > result;
    for ( int i = 0; i <= order; i++ )
      for ( auto rest : exponents( order - i, dim - 1 ) ) {
        rest.push_back( i );
        result.push_back( rest );
      }
    return result;
  }

  double monomial( const std::vector< int >& a, const Eigen::VectorXd& x )
  {
    double m = 1;
    for ( size_t d = 0; d < a.size(); d++ )
      m *= std::pow( x[d], a[d] );
    return m;
  }

  double monomialDerivative( const std::vector< int >& a, const Eigen::VectorXd& x, int i )
  {
    if ( a[i] == 0 )
      return 0;
    double m = a[i] * std::pow( x[i], a[i] - 1 );
    for ( size_t d = 0; d < a.size(); d++ )
      if ( int( d ) != i )
        m *= std::pow( x[d], a[d] );
    return m;
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

  Eigen::VectorXd point( int dim, double offset )
  {
    Eigen::VectorXd x( dim );
    for ( int d = 0; d < dim; d++ )
      x[d] = offset + 0.13 * ( d + 1 );
    return x;
  }

  // a regular grid of n^dim kernels with spacing 1
  template < typename Kernel >
  struct KernelGrid {
    std::vector< Eigen::VectorXd >                                          centers;
    std::vector< std::unique_ptr< Kernel > >                                kernels;
    std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* > pointers;

    KernelGrid( int dim, int n, double supportRadius )
    {
      const int nTotal = static_cast< int >( std::pow( n, dim ) );
      centers.reserve( nTotal ); // the kernels keep pointers to their centres
      for ( int k = 0; k < nTotal; k++ ) {
        Eigen::VectorXd c( dim );
        for ( int d = 0, r = k; d < dim; d++, r /= n )
          c[d] = r % n;
        centers.push_back( c );
        kernels.emplace_back( std::make_unique< Kernel >( centers.back().data(), dim, supportRadius ) );
        pointers.push_back( kernels.back().get() );
      }
    }
  };

  // The visitors of a particle's stiffness, fed by hand with a node description of 5 dofs per node in 3D (no particle
  // of Marmot has that many, so it is not a family the visitors know): the assembled stiffness is the chain rule of
  // StiffnessVisitor, block by block, and the bound of |K| W is never below it. A particle with more dofs per node than
  // the visitors are compiled for is refused.
  void testStiffnessVisitorsOnAnyDofsPerNode()
  {
    using namespace Marmot::Discretizations::Meshfree;
    constexpr int                    nDim = 3, n = 5, nNodes = 4, nDof = n * nNodes;
    std::mt19937                     gen( 7 );
    std::uniform_real_distribution<> u( -1.0, 1.0 );
    const auto                       random = [&]( int rows, int cols ) {
      Eigen::MatrixXd m( rows, cols );
      for ( int i = 0; i < m.size(); i++ )
        m( i ) = u( gen );
      return m;
    };
    const Eigen::MatrixXd N = random( 1, nNodes ), dN_dY = random( nDim, nNodes );
    Eigen::MatrixXd       d2N_dYdY( nDim * nDim, nNodes );
    for ( int B = 0; B < nNodes; B++ )
      for ( int l = 0; l < nDim; l++ )
        for ( int m = l; m < nDim; m++ )
          d2N_dYdY( l * nDim + m, B ) = d2N_dYdY( m * nDim + l, B ) = u( gen );
    std::vector< Eigen::MatrixXd > a, b, h;
    for ( int A = 0; A < nNodes; A++ ) {
      a.push_back( random( n, n ) );
      b.push_back( random( n, nDim * n ) );
      h.push_back( random( n, nDim * nDim * n ) );
    }

    Eigen::MatrixXd                K = Eigen::MatrixXd::Zero( nDof, nDof ), bound = Eigen::MatrixXd::Zero( nDof, 1 );
    const Eigen::MatrixXd          W = Eigen::MatrixXd::Ones( nDof, 1 );
    StiffnessAssembler< nDim >     assembler( K.data() );
    AbsoluteStiffnessBound< nDim > absoluteStiffness( W.data(), 1, bound.data() );
    for ( StiffnessVisitor* visitor : std::vector< StiffnessVisitor* >{ &assembler, &absoluteStiffness } ) {
      visitor->visitShapeFunctions( nDim, n, nNodes, N.data(), dN_dY.data(), d2N_dYdY.data() );
      for ( int A = 0; A < nNodes; A++ )
        visitor->visitStiffnessOfNode( A, a[A].data(), b[A].data(), h[A].data() );
    }

    Eigen::MatrixXd byHand = Eigen::MatrixXd::Zero( nDof, nDof );
    for ( int A = 0; A < nNodes; A++ )
      for ( int B = 0; B < nNodes; B++ ) {
        Eigen::MatrixXd K_AB = a[A] * N( 0, B );
        for ( int l = 0; l < nDim; l++ )
          K_AB += b[A].middleCols( n * l, n ) * dN_dY( l, B );
        for ( int lm = 0; lm < nDim * nDim; lm++ )
          K_AB += h[A].middleCols( n * lm, n ) * d2N_dYdY( lm, B );
        byHand.block( n * A, n * B, n, n ) = K_AB;
      }
    throwExceptionOnFailure( ( K - byHand ).cwiseAbs().maxCoeff() <= 1e-14 * byHand.cwiseAbs().maxCoeff(),
                             "StiffnessAssembler: the stiffness of 5 dofs per node is not the chain rule" );
    const Eigen::MatrixXd exact = byHand.cwiseAbs() * W;
    throwExceptionOnFailure( ( bound - exact ).minCoeff() >= -1e-12 * exact.maxCoeff(),
                             "AbsoluteStiffnessBound: the bound of 5 dofs per node lies below |K| W" );

    for ( StiffnessVisitor* visitor : std::vector< StiffnessVisitor* >{ &assembler, &absoluteStiffness } ) {
      bool refused = false;
      try {
        visitor->visitShapeFunctions( nDim, maxDofsPerNode + 1, nNodes, N.data(), dN_dY.data(), nullptr );
      }
      catch ( const std::invalid_argument& ) {
        refused = true;
      }
      throwExceptionOnFailure( refused, "a stiffness visitor reads a node with more dofs than it is compiled for" );
    }
  }

} // namespace

void testMonomialBasis()
{
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 0; order <= 3; order++ ) {
      const auto a = exponents( order, dim );
      const int  n = Math::computeSizeOfMonomialBasisVector( order, dim );
      throwExceptionOnFailure( n == int( a.size() ),
                               Marmot::Utilities::MakeString()
                                 << "size of the basis of order " << order << " in " << dim << "D" );

      const Eigen::VectorXd x = point( dim, 0.7 );
      Eigen::VectorXd       H( n );
      Eigen::MatrixXd       dH( n, dim );
      Math::computeMonomialBasis( order, x, H );
      Math::computeMonomialBasisGradient( order, x, dH );

      // the basis is ordered as the exponent tuples above; compare with the direct evaluation
      for ( int k = 0; k < n; k++ ) {
        throwExceptionOnFailure( checkIfEqual( H[k], monomial( a[k], x ), 1e-14 ),
                                 Marmot::Utilities::MakeString()
                                   << "monomial " << k << " of order " << order << " in " << dim << "D" );
        for ( int i = 0; i < dim; i++ )
          throwExceptionOnFailure( checkIfEqual( dH( k, i ), monomialDerivative( a[k], x, i ), 1e-13 ),
                                   Marmot::Utilities::MakeString()
                                     << "derivative " << i << " of monomial " << k << " of order " << order << " in "
                                     << dim << "D: " << dH( k, i ) << " vs " << monomialDerivative( a[k], x, i ) );
      }
    }
}

template < typename Kernel >
void checkKernel( const std::string& name )
{
  for ( int dim = 1; dim <= 3; dim++ ) {
    Eigen::VectorXd c = point( dim, 0.2 );
    const double    r = 1.3;
    Kernel          kernel( c.data(), dim, r );

    throwExceptionOnFailure( kernel.getCenterCoordinates() == c.data(), name + ": center coordinates" );
    throwExceptionOnFailure( kernel.computeKernelFunction( c.data() ) > 0, name + ": positive at the center" );

    Eigen::VectorXd lower( dim ), upper( dim );
    kernel.getBoundingBox( lower.data(), upper.data() );
    throwExceptionOnFailure( ( lower - ( c.array() - r ).matrix() ).norm() < 1e-14 &&
                               ( upper - ( c.array() + r ).matrix() ).norm() < 1e-14,
                             name + ": bounding box" );

    // outside the box: no support, zero value
    Eigen::VectorXd outside = c;
    outside[dim - 1] += 1.01 * r;
    throwExceptionOnFailure( !kernel.isInSupport( outside.data() ) &&
                               kernel.computeKernelFunction( outside.data() ) == 0,
                             name + ": support" );

    // the gradient is the derivative of the value
    for ( double s : { 0.1, 0.45, 0.8 } ) {
      Eigen::VectorXd x = c;
      for ( int d = 0; d < dim; d++ )
        x[d] += s * r * ( d % 2 ? -1 : 1 ) / ( d + 1 );
      throwExceptionOnFailure( kernel.isInSupport( x.data() ), name + ": support inside" );

      Eigen::VectorXd grad( dim );
      kernel.computeKernelFunctionGradient( x.data(), grad.data() );
      for ( int i = 0; i < dim; i++ ) {
        const double    h  = 1e-6;
        Eigen::VectorXd xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        const double fd = ( kernel.computeKernelFunction( xp.data() ) - kernel.computeKernelFunction( xm.data() ) ) /
                          ( 2 * h );
        throwExceptionOnFailure( std::abs( grad[i] - fd ) < 1e-7 * ( 1 + std::abs( fd ) ),
                                 Marmot::Utilities::MakeString()
                                   << name << ": gradient " << i << " in " << dim << "D" );
      }
    }

    // moving the kernel moves its support
    Eigen::VectorXd target = c;
    target[0] += 5.0;
    kernel.moveTo( target.data() );
    throwExceptionOnFailure( kernel.isInSupport( target.data() ) && !kernel.isInSupport( point( dim, 0.2 ).data() ),
                             name + ": moveTo" );
  }
}

void testKernelFunctions()
{
  checkKernel< Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed >( "BSpline2ndOrderBoxed" );
  checkKernel< Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed >( "BSpline3rdOrderBoxed" );
}

// reproduction of the complete polynomial basis of the requested order by the values, and of its derivatives by
// the gradients; for the explicit (direct) derivatives, also consistency with the shape function values
template < typename Kernel, typename Approximation >
void checkApproximation( const std::string& name, int dim, int order, double supportRadius, bool directDerivatives )
{
  KernelGrid< Kernel > grid( dim, 6, supportRadius );
  const Approximation  approximation( dim, order );
  const int            nNodes = static_cast< int >( grid.pointers.size() );
  const auto           a      = exponents( order, dim );

  auto values = [&]( const Eigen::VectorXd& x ) {
    Eigen::VectorXd N( nNodes );
    approximation.computeShapeFunctions( x.data(), grid.pointers, N.data() );
    return N;
  };

  for ( double offset : { 2.1, 2.55, 3.0 } ) {
    const Eigen::VectorXd x = point( dim, offset );
    Eigen::VectorXd       N( nNodes );
    Eigen::MatrixXd       dN( dim, nNodes );
    approximation.computeShapeFunctionsAndGradients( x.data(), grid.pointers, N.data(), dN.data() );

    throwExceptionOnFailure( ( N - values( x ) ).norm() < 1e-12, name + ": values of both evaluations differ" );

    for ( const auto& alpha : a ) {
      double          reproduced         = 0;
      Eigen::VectorXd reproducedGradient = Eigen::VectorXd::Zero( dim );
      for ( int A = 0; A < nNodes; A++ ) {
        reproduced += N[A] * monomial( alpha, grid.centers[A] );
        reproducedGradient += dN.col( A ) * monomial( alpha, grid.centers[A] );
      }
      throwExceptionOnFailure( std::abs( reproduced - monomial( alpha, x ) ) < 1e-10,
                               Marmot::Utilities::MakeString()
                                 << name << ", order " << order << " in " << dim << "D: monomial not reproduced" );
      for ( int i = 0; i < dim; i++ )
        throwExceptionOnFailure( std::abs( reproducedGradient[i] - monomialDerivative( alpha, x, i ) ) < 1e-9,
                                 Marmot::Utilities::MakeString()
                                   << name << ", order " << order << " in " << dim
                                   << "D: derivative of a monomial not reproduced, " << reproducedGradient[i] << " vs "
                                   << monomialDerivative( alpha, x, i ) );
    }

    if ( directDerivatives )
      for ( int i = 0; i < dim; i++ ) {
        const double    h  = 1e-6;
        Eigen::VectorXd xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        const Eigen::VectorXd fd = ( values( xp ) - values( xm ) ) / ( 2 * h );
        throwExceptionOnFailure( ( dN.row( i ).transpose() - fd ).norm() < 1e-6 * ( 1 + fd.norm() ),
                                 Marmot::Utilities::MakeString()
                                   << name << ", order " << order << " in " << dim
                                   << "D: gradients are not the derivatives of the values" );
      }
  }
}

void testReproducingKernelApproximation()
{
  using K2 = Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed;
  using K3 = Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed;
  using RK = Marmot::Discretizations::Meshfree::ReproducingKernelApproximation;
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      const double r = order == 1 ? 1.6 : 2.4;
      checkApproximation< K2, RK >( "RK, 2nd order B-spline", dim, order, r, true );
      checkApproximation< K3, RK >( "RK, 3rd order B-spline", dim, order, r, true );
    }
}

void testImplicitGradientReproducingKernelApproximation()
{
  using K2  = Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed;
  using K3  = Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed;
  using IRK = Marmot::Discretizations::Meshfree::ReproducingKernelApproximationImplicit;
  // implicit gradients reproduce the derivatives of the basis, but are not the derivatives of the values
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      const double r = order == 1 ? 1.6 : 2.4;
      checkApproximation< K2, IRK >( "implicit RK, 2nd order B-spline", dim, order, r, false );
      checkApproximation< K3, IRK >( "implicit RK, 3rd order B-spline", dim, order, r, false );
    }
}

// the moment matrix assembly is protected: expose it for the test, with the kernel evaluation it takes its input from
struct MomentMatrixAccess : Marmot::Discretizations::Meshfree::ReproducingKernelApproximation {
  using ReproducingKernelApproximation::conditionOfMomentMatrix;
  using ReproducingKernelApproximation::ReproducingKernelApproximation;

  // the moment matrix at x, and its gradients if dM_dx is given
  void momentMatrix( const Eigen::VectorXd&                                                         x,
                     const std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* >& kernels,
                     int                                                                            order,
                     MomentMatrix&                                                                  M,
                     MomentMatrixGradient*                                                          dM_dx ) const
  {
    std::vector< double > phi( kernels.size() ), dPhi_dx( _dim * kernels.size() );
    evaluateKernelFunctions( x.data(), kernels, phi.data(), dPhi_dx.data() );
    assembleMomentMatrix( x.data(), kernels, phi.data(), dPhi_dx.data(), order, M, dM_dx );
  }
};

void testMomentMatrixGradient()
{
  using RK = MomentMatrixAccess;
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      KernelGrid< Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed > grid( dim, 6, 2.4 );
      const RK                                                                            rk( dim, order );
      const Eigen::VectorXd                                                               x = point( dim, 2.3 );
      RK::MomentMatrix                                                                    M, MWithoutGradients;
      RK::MomentMatrixGradient                                                            dM;
      rk.momentMatrix( x, grid.pointers, order, M, &dM );
      rk.momentMatrix( x, grid.pointers, order, MWithoutGradients, nullptr );
      throwExceptionOnFailure( ( M - MWithoutGradients ).norm() < 1e-14 * M.norm(),
                               "moment matrix with and without its gradients differs" );
      for ( int i = 0; i < dim; i++ ) {
        const double    h  = 1e-6;
        Eigen::VectorXd xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        RK::MomentMatrix Mp, Mm;
        rk.momentMatrix( xp, grid.pointers, order, Mp, nullptr );
        rk.momentMatrix( xm, grid.pointers, order, Mm, nullptr );
        const Eigen::MatrixXd fd = ( Mp - Mm ) / ( 2 * h );
        throwExceptionOnFailure( ( dM[i] - fd ).norm() < 1e-6 * ( 1 + fd.norm() ),
                                 Marmot::Utilities::MakeString()
                                   << "gradient " << i << " of the moment matrix of order " << order << " in " << dim
                                   << "D" );
      }
    }
}

void testCompletenessOrderIsReducedForFewNodes()
{
  // two kernels in 1D cannot carry a quadratic basis: the order drops to 1
  KernelGrid< Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed > grid( 1, 2, 1.6 );
  const Marmot::Discretizations::Meshfree::ReproducingKernelApproximation             approximation( 1, 2 );
  const double                                                                        x = 0.4;
  double                                                                              N[2], dN[2];
  approximation.computeShapeFunctionsAndGradients( &x, grid.pointers, N, dN );
  throwExceptionOnFailure( checkIfEqual( N[0] + N[1], 1.0, 1e-14 ), "partition of unity with a reduced order" );
  throwExceptionOnFailure( std::abs( N[0] * 0 + N[1] * 1 - x ) < 1e-12,
                           "two nodes in 1D still reproduce a linear field" );
}

// the evaluation is of a fixed largest size: a basis beyond it, a dimension beyond 3 or a negative order is refused
void testUnsupportedBasisIsRefused()
{
  using RK = Marmot::Discretizations::Meshfree::ReproducingKernelApproximation;
  throwExceptionOnFailure( !throws( [] { RK( 3, 3 ); } ), "the cubic basis in 3D is supported" );
  for ( const auto& [dim, order] : { std::pair( 3, 4 ), std::pair( 4, 1 ), std::pair( 2, -1 ) } )
    throwExceptionOnFailure( throws( [dim = dim, order = order] { RK( dim, order ); } ),
                             Marmot::Utilities::MakeString()
                               << "dimension " << dim << " and order " << order << " must be refused" );
}

// a singular moment matrix must throw instead of silently losing the partition of unity: three collinear nodes in 2D
// pass the node count check of a linear basis, and a point outside all supports has no kernel at all
void testSingularMomentMatrixThrows()
{
  using K = Marmot::Discretizations::Meshfree::KernelFunctionBSpline2ndOrderBoxed;
  std::vector< Eigen::VectorXd >                                          centers = { Eigen::Vector2d( 0, 0 ),
                                                                                      Eigen::Vector2d( 1, 0 ),
                                                                                      Eigen::Vector2d( 2, 0 ) };
  std::vector< std::unique_ptr< K > >                                     kernels;
  std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* > pointers;
  for ( auto& c : centers ) {
    kernels.emplace_back( std::make_unique< K >( c.data(), 2, 1.6 ) );
    pointers.push_back( kernels.back().get() );
  }
  const Marmot::Discretizations::Meshfree::ReproducingKernelApproximation         rk( 2, 1 );
  const Marmot::Discretizations::Meshfree::ReproducingKernelApproximationImplicit irk( 2, 1 );

  for ( const Eigen::Vector2d& x : { Eigen::Vector2d( 1.1, 0.2 ), Eigen::Vector2d( 10, 10 ) } ) {
    double N[3], dN[6];
    throwExceptionOnFailure( throws( [&]() { rk.computeShapeFunctions( x.data(), pointers, N ); } ) &&
                               throws( [&]() { rk.computeShapeFunctionsAndGradients( x.data(), pointers, N, dN ); } ) &&
                               throws( [&]() { irk.computeShapeFunctionsAndGradients( x.data(), pointers, N, dN ); } ),
                             Marmot::Utilities::MakeString()
                               << "a singular moment matrix at " << x.transpose() << " must throw" );
  }
}

// One layer of kernels on a lattice of spacing 1, support radius 1.6, every other centre lifted by a small jitter out
// of the plane z = 0 (a wall of particles), and a point at z = 0.8 above it (in the void next to the wall): the
// covering kernels lie nearly in one plane, on one side of the point. Reproducing the linear field z there needs shape
// functions of the size of the distance over the jitter: the moment matrix of the linear basis is nearly singular, and
// singular without a jitter.
struct OneSidedKernels {
  using K = Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed;
  std::vector< Eigen::VectorXd >                                          centers;
  std::vector< std::unique_ptr< K > >                                     kernels;
  std::vector< const Marmot::Discretizations::Meshfree::KernelFunction* > pointers;
  Eigen::Vector3d                                                         point;

  explicit OneSidedKernels( double jitter, double length = 1.0 ) : point( length * Eigen::Vector3d( 2.1, 1.8, 0.8 ) )
  {
    centers.reserve( 25 ); // the kernels keep pointers to their centres
    for ( int j = 0; j < 5; j++ )
      for ( int i = 0; i < 5; i++ ) {
        centers.push_back( length * Eigen::Vector3d( i, j, ( i + j ) % 2 ? jitter : 0.0 ) );
        kernels.emplace_back( std::make_unique< K >( centers.back().data(), 3, length * 1.6 ) );
        pointers.push_back( kernels.back().get() );
      }
  }
};

// The optional limit on the condition of the moment matrix: at the one-sided point, the linear shape functions of the
// unlimited approximation are large; above the limit the point gets the Shepard functions (order 0): a partition of
// unity of values in [0, 1], with gradients of the size of the inverse support radius. The condition is that of the
// moment matrix scaled to a unit diagonal, the same in any unit of length.
void testConditionLimitOfTheMomentMatrix()
{
  using RK  = Marmot::Discretizations::Meshfree::ReproducingKernelApproximation;
  using IRK = Marmot::Discretizations::Meshfree::ReproducingKernelApproximationImplicit;
  const OneSidedKernels oneSided( 0.02 );
  const int             n = static_cast< int >( oneSided.pointers.size() );

  const MomentMatrixAccess access( 3, 1 );
  RK::MomentMatrix         M;
  access.momentMatrix( oneSided.point, oneSided.pointers, 1, M, nullptr );
  const double condition = access.conditionOfMomentMatrix( M );
  throwExceptionOnFailure( condition > 1e4,
                           Marmot::Utilities::MakeString()
                             << "the one-sided point must be nearly singular, condition " << condition );

  const OneSidedKernels inMillimetres( 0.02, 1000.0 );
  RK::MomentMatrix      MInMillimetres;
  access.momentMatrix( inMillimetres.point, inMillimetres.pointers, 1, MInMillimetres, nullptr );
  throwExceptionOnFailure( checkIfEqual( access.conditionOfMomentMatrix( MInMillimetres ),
                                         condition,
                                         1e-6 * condition ),
                           "the condition of the scaled moment matrix depends on the unit of length" );

  Eigen::VectorXd N( n ), NLimited( n );
  Eigen::MatrixXd dN( 3, n ), dNLimited( 3, n );
  RK( 3, 1 ).computeShapeFunctionsAndGradients( oneSided.point.data(), oneSided.pointers, N.data(), dN.data() );
  RK( 3, 1, 1e4 )
    .computeShapeFunctionsAndGradients( oneSided.point.data(), oneSided.pointers, NLimited.data(), dNLimited.data() );
  throwExceptionOnFailure( N.cwiseAbs().maxCoeff() > 10.0,
                           Marmot::Utilities::MakeString()
                             << "without a limit, the one-sided point must have large shape functions, max |N| "
                             << N.cwiseAbs().maxCoeff() );

  double sumOfShepardKernels = 0;
  for ( int A = 0; A < n; A++ )
    sumOfShepardKernels += oneSided.pointers[A]->computeKernelFunction( oneSided.point.data() );
  for ( int A = 0; A < n; A++ ) {
    const double phi_A = oneSided.pointers[A]->computeKernelFunction( oneSided.point.data() );
    throwExceptionOnFailure( std::abs( NLimited( A ) - phi_A / sumOfShepardKernels ) < 1e-14,
                             "above the limit, the shape functions are the Shepard functions" );
  }
  throwExceptionOnFailure( checkIfEqual( NLimited.sum(), 1.0, 1e-14 ) && NLimited.minCoeff() >= 0.0 &&
                             dNLimited.cwiseAbs().maxCoeff() * 1.6 < 10.0 && dNLimited.rowwise().sum().norm() < 1e-12,
                           Marmot::Utilities::MakeString()
                             << "above the limit: a partition of unity in [0, 1] with bounded gradients, max |dN| "
                             << dNLimited.cwiseAbs().maxCoeff() );
  for ( int i = 0; i < 3; i++ ) {
    const double    h  = 1e-6;
    Eigen::Vector3d xp = oneSided.point, xm = oneSided.point;
    xp[i] += h;
    xm[i] -= h;
    Eigen::VectorXd Np( n ), Nm( n );
    RK( 3, 1, 1e4 ).computeShapeFunctions( xp.data(), oneSided.pointers, Np.data() );
    RK( 3, 1, 1e4 ).computeShapeFunctions( xm.data(), oneSided.pointers, Nm.data() );
    throwExceptionOnFailure( ( dNLimited.row( i ).transpose() - ( Np - Nm ) / ( 2 * h ) ).norm() < 1e-6,
                             "above the limit, the gradients are those of the Shepard functions" );
  }
  // the implicit approximation chooses the order as the explicit one, and at the order 0 its gradients are the
  // Shepard gradients as well (its implicit gradients of order 0 would vanish identically)
  Eigen::VectorXd NImplicit( n );
  Eigen::MatrixXd dNImplicit( 3, n );
  IRK( 3, 1, 1e4 )
    .computeShapeFunctionsAndGradients( oneSided.point.data(), oneSided.pointers, NImplicit.data(), dNImplicit.data() );
  throwExceptionOnFailure( ( NImplicit - NLimited ).norm() < 1e-14,
                           "the implicit approximation chooses the order as the explicit one" );
  throwExceptionOnFailure( dNLimited.cwiseAbs().maxCoeff() > 0.1 && ( dNImplicit - dNLimited ).norm() < 1e-12,
                           Marmot::Utilities::MakeString()
                             << "above the limit, the implicit gradients must be the Shepard gradients, max |dN| "
                             << dNImplicit.cwiseAbs().maxCoeff() << ", difference "
                             << ( dNImplicit - dNLimited ).norm() );

  // a healthy point is evaluated bit for bit as without a limit, by both approximations
  KernelGrid< Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed > grid( 3, 6, 1.6 );
  const int             nGrid = static_cast< int >( grid.pointers.size() );
  const Eigen::VectorXd x     = point( 3, 2.55 );
  for ( bool implicit : { false, true } ) {
    std::unique_ptr< RK > unlimited = implicit ? std::make_unique< IRK >( 3, 1 ) : std::make_unique< RK >( 3, 1 );
    std::unique_ptr< RK > limited   = implicit ? std::make_unique< IRK >( 3, 1, 1e4 )
                                               : std::make_unique< RK >( 3, 1, 1e4 );
    Eigen::VectorXd       NGrid( nGrid ), NGridLimited( nGrid );
    Eigen::MatrixXd       dNGrid( 3, nGrid ), dNGridLimited( 3, nGrid );
    unlimited->computeShapeFunctionsAndGradients( x.data(), grid.pointers, NGrid.data(), dNGrid.data() );
    limited->computeShapeFunctionsAndGradients( x.data(), grid.pointers, NGridLimited.data(), dNGridLimited.data() );
    throwExceptionOnFailure( NGrid == NGridLimited && dNGrid == dNGridLimited,
                             "a healthy point is not evaluated bit for bit as without a limit" );
  }

  // a singular moment matrix of the linear basis (all covering kernels in one plane) gets the order 0 with a limit
  const OneSidedKernels coplanar( 0.0 );
  throwExceptionOnFailure( throws( [&] {
                             RK( 3, 1 ).computeShapeFunctions( coplanar.point.data(), coplanar.pointers, N.data() );
                           } ) &&
                             !throws( [&] {
                               RK( 3, 1, 1e4 )
                                 .computeShapeFunctions( coplanar.point.data(), coplanar.pointers, N.data() );
                             } ) &&
                             checkIfEqual( N.sum(), 1.0, 1e-14 ),
                           "a singular moment matrix of the linear basis must throw without a limit and get the "
                           "Shepard functions with one" );

  // an exactly singular matrix is infinitely ill-conditioned, also when its factorization ends on a zero pivot
  for ( const Eigen::Matrix2d& singular :
        { Eigen::Matrix2d{ { 1.0, 1.0 }, { 1.0, 1.0 } }, Eigen::Matrix2d{ { 4.0, 2.0 }, { 2.0, 1.0 } } } )
    throwExceptionOnFailure( std::isinf( MomentMatrixAccess::conditionOfMomentMatrix( singular ) ),
                             Marmot::Utilities::MakeString()
                               << "a singular matrix must have an infinite condition, not "
                               << MomentMatrixAccess::conditionOfMomentMatrix( singular ) );

  for ( double limit : { 0.5, -1.0, std::numeric_limits< double >::quiet_NaN() } )
    throwExceptionOnFailure( throws( [limit] { RK( 3, 1, limit ); } ) && throws( [limit] { IRK( 3, 1, limit ); } ),
                             Marmot::Utilities::MakeString() << "the condition limit " << limit << " must be refused" );
}

// The condition of the scaled moment matrix of the linear basis on a regular lattice of kernels (spacing 1): inside,
// and at a face, an edge and a corner of a body of particles whose kernels sit at the particle centres, half a spacing
// beyond the outermost centres. These are the bounds the documentation quotes for the choice of the limit.
void testConditionOfTheMomentMatrixOnALattice()
{
  using K = Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed;
  const MomentMatrixAccess access( 3, 1 );
  for ( const auto& [supportRadius, bound] : { std::pair( 1.6, 1e4 ), std::pair( 2.4, 1e2 ) } ) {
    KernelGrid< K > grid( 3, 6, supportRadius );
    for ( const Eigen::Vector3d& x : { Eigen::Vector3d( 2.55, 2.45, 2.6 ),
                                       Eigen::Vector3d( -0.5, 2.45, 2.6 ),
                                       Eigen::Vector3d( -0.5, -0.5, 2.6 ),
                                       Eigen::Vector3d( -0.5, -0.5, -0.5 ) } ) {
      MomentMatrixAccess::MomentMatrix M;
      access.momentMatrix( x, grid.pointers, 1, M, nullptr );
      const double condition = MomentMatrixAccess::conditionOfMomentMatrix( M );
      throwExceptionOnFailure( condition < bound,
                               Marmot::Utilities::MakeString()
                                 << "support radius " << supportRadius << ", point " << x.transpose()
                                 << ": the condition " << condition << " exceeds " << bound );
    }
  }
}

// With a limit, the order of a point is lowered one order at a time: near the corner of a body, the quadratic basis
// exceeds the limit where the linear one does not, and the point keeps the linear completeness; a point where the
// quadratic basis is below the limit keeps it.
void testOrderIsLoweredOneAtATime()
{
  using K  = Marmot::Discretizations::Meshfree::KernelFunctionBSpline3rdOrderBoxed;
  using RK = Marmot::Discretizations::Meshfree::ReproducingKernelApproximation;
  KernelGrid< K >          grid( 3, 6, 2.4 );
  const int                n = static_cast< int >( grid.pointers.size() );
  const MomentMatrixAccess access( 3, 2 );
  const double             limit = 1e3;

  const Eigen::Vector3d            nearTheCorner( -0.25, -0.25, -0.25 );
  MomentMatrixAccess::MomentMatrix M1, M2;
  access.momentMatrix( nearTheCorner, grid.pointers, 1, M1, nullptr );
  access.momentMatrix( nearTheCorner, grid.pointers, 2, M2, nullptr );
  throwExceptionOnFailure( MomentMatrixAccess::conditionOfMomentMatrix( M2 ) > limit &&
                             MomentMatrixAccess::conditionOfMomentMatrix( M1 ) < limit,
                           "near the corner, the quadratic basis must exceed the limit and the linear one not" );

  Eigen::VectorXd N( n ), NLinear( n );
  Eigen::MatrixXd dN( 3, n ), dNLinear( 3, n );
  RK( 3, 2, limit ).computeShapeFunctionsAndGradients( nearTheCorner.data(), grid.pointers, N.data(), dN.data() );
  RK( 3, 1 ).computeShapeFunctionsAndGradients( nearTheCorner.data(), grid.pointers, NLinear.data(), dNLinear.data() );
  Eigen::Vector3d reproduced = Eigen::Vector3d::Zero();
  for ( int A = 0; A < n; A++ )
    reproduced += N( A ) * grid.centers[A];
  throwExceptionOnFailure( ( N - NLinear ).norm() < 1e-14 && ( dN - dNLinear ).norm() < 1e-12 &&
                             ( reproduced - nearTheCorner ).norm() < 1e-12,
                           "above the limit of the quadratic basis, the point must get the linear one" );

  const Eigen::Vector3d atTheFace( 0.0, 2.45, 2.6 );
  Eigen::VectorXd       NQuadratic( n );
  RK( 3, 2, limit ).computeShapeFunctions( atTheFace.data(), grid.pointers, N.data() );
  RK( 3, 2 ).computeShapeFunctions( atTheFace.data(), grid.pointers, NQuadratic.data() );
  throwExceptionOnFailure( N == NQuadratic, "below the limit, the point must keep the quadratic basis" );
}

// the inverse isoparametric map and the point location test of a distorted cell (Newton's method), and of a box cell
template < int nDim, int nNodes >
void checkLagrangianCellInverseMap( const std::vector< double >& nodes, bool box )
{
  using Geometry = Marmot::Discretizations::Cells::LagrangianCellGeometry< nDim, nNodes >;
  using Xi       = typename Geometry::XiSized;
  const Geometry geometry( nodes.data() );
  const auto     X = Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > >( nodes.data() );

  for ( const double a : { -0.9, -0.3, 0.0, 0.45, 0.8 } ) {
    Xi xi;
    for ( int i = 0; i < nDim; i++ )
      xi( i ) = a + 0.07 * i;
    const Xi x = X * geometry.N( xi ).transpose();
    throwExceptionOnFailure( ( geometry.findReferenceCoordinate( x ) - xi ).norm() < 1e-10,
                             Marmot::Utilities::MakeString() << nDim << "D cell: inverse map of " << xi.transpose() );
    throwExceptionOnFailure( geometry.isCoordinateInCell( x.data() ),
                             Marmot::Utilities::MakeString()
                               << nDim << "D cell: a point at " << xi.transpose() << " is in the cell" );
  }

  // a point of the bounding box outside a distorted cell: xi_0 = 1.2 (inside its bounding box for these cells)
  Xi xi      = Xi::Zero();
  xi( 0 )    = 1.2;
  const Xi x = X * geometry.N( xi ).transpose();
  double   lo[nDim], hi[nDim];
  geometry.getBoundingBox( lo, hi );
  bool inBox = true;
  for ( int i = 0; i < nDim; i++ )
    inBox = inBox && x( i ) >= lo[i] && x( i ) < hi[i];
  if ( !box && inBox )
    throwExceptionOnFailure( !geometry.isCoordinateInCell( x.data() ),
                             Marmot::Utilities::MakeString()
                               << nDim << "D cell: a point outside the distorted cell is not in it" );
}

void testLagrangianCellInverseMap()
{
  checkLagrangianCellInverseMap< 2, 4 >( { 0, 0, 2, 0, 2, 1, 0, 1 }, true );
  checkLagrangianCellInverseMap< 2, 4 >( { 0, 0, 2, 0.3, 2.6, 1.4, -0.2, 1 }, false );
  checkLagrangianCellInverseMap< 3, 8 >( { 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 1, 0, 0, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1, 1 },
                                         true );
  checkLagrangianCellInverseMap< 3, 8 >( { 0, 0, 0, 2,   0.2, 0,   2.4, 1.3, 0.1, -0.1, 1, 0,
                                           0, 0, 1, 2.2, 0,   1.2, 2.5, 1.2, 1.3, 0,    1, 1.1 },
                                         false );
}

// the subdomains of a particle domain are deformed about the centroid of the parent, so that they tile the
// deformed parent (about their own centroids, they would overlap or leave gaps for any F != I)
void testSubdomainsTileTheDeformedDomain()
{
  using Domain               = ParticleDomain< 2, 4 >;
  const double    vertices[] = { -1, -1, 1, -1, 1, 1, -1, 1 };
  Domain          main( vertices, 8, Domain::DeformationGradient );
  auto            subdomains = main.uniformSubdivided();
  Eigen::Matrix2d F;
  F << 2, 0.3, 0.1, 1;
  const Eigen::Vector2d u( 0.4, -0.2 );

  main.acceptStateAndPosition( F, u );
  double volume = 0;
  for ( auto& sd : subdomains ) {
    const Eigen::Matrix< double, 2, 4 > X = sd.getGeometryDeformedVertexCoordinates(); // undeformed, not accepted
    sd.acceptStateAndPosition( F, u );
    const Eigen::Matrix< double, 2, 4 > x = sd.getGeometryDeformedVertexCoordinates();
    throwExceptionOnFailure( ( x - ( ( F * X ).colwise() + u ) ).norm() < 1e-13,
                             "a subdomain follows the homogeneous map of the parent about the parent centroid" );
    volume += Domain::getVolumeFromVertices( x );
  }
  throwExceptionOnFailure( checkIfEqual( volume,
                                         Domain::getVolumeFromVertices( main.getGeometryDeformedVertexCoordinates() ),
                                         1e-12 ),
                           "the deformed subdomains tile the deformed domain" );
}

// the factories: a Lagrangian and a B-spline cell may share a name (the B-spline registration checked the wrong map
// for duplicates), names are case-insensitive, and an unknown name is reported by name
void testFactories()
{

  const std::string name = "Test/FactoryCell";
  throwExceptionOnFailure( Marmot::Registry::CellFactory::registerCell( name,
                                                                        []( int, const double*, int )
                                                                          -> Marmot::Discretizations::Cells::Cell* {
                                                                          return nullptr;
                                                                        } ),
                           "Lagrangian cell registration" );
  throwExceptionOnFailure( Marmot::Registry::CellFactory::
                             registerBSplineCell( name,
                                                  []( int, const double*, int, const double*, int )
                                                    -> Marmot::Discretizations::Cells::Cell* { return nullptr; } ),
                           "a B-spline cell may have the name of a Lagrangian cell" );
  throwExceptionOnFailure( Marmot::Registry::CellFactory::createCell( "TEST/FACTORYCELL", 1, nullptr, 0 ) == nullptr &&
                             Marmot::Registry::CellFactory::createBSplineCell( "test/factorycell",
                                                                               1,
                                                                               nullptr,
                                                                               0,
                                                                               nullptr,
                                                                               0 ) == nullptr,
                           "cells are created by their case-insensitive name" );

  const std::string unknown  = "No/Such/CellElement";
  bool              reported = false;
  try {
    Marmot::Registry::CellElementFactory::createCellElement( unknown, 42, nullptr, 0, "gauss", 1 );
  }
  catch ( const std::invalid_argument& e ) {
    reported = std::string( e.what() ).find( unknown ) != std::string::npos;
  }
  throwExceptionOnFailure( reported, "an unknown cell element is reported by its name" );

  bool threw = false;
  try {
    Marmot::Registry::MaterialPointFactory::createMaterialPoint( "No/Such/MaterialPoint", 1, nullptr, 0, 1.0 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "an unknown material point must throw" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testMonomialBasis,
                                                               testKernelFunctions,
                                                               testReproducingKernelApproximation,
                                                               testImplicitGradientReproducingKernelApproximation,
                                                               testMomentMatrixGradient,
                                                               testCompletenessOrderIsReducedForFewNodes,
                                                               testUnsupportedBasisIsRefused,
                                                               testSingularMomentMatrixThrows,
                                                               testConditionLimitOfTheMomentMatrix,
                                                               testConditionOfTheMomentMatrixOnALattice,
                                                               testOrderIsLoweredOneAtATime,
                                                               testLagrangianCellInverseMap,
                                                               testSubdomainsTileTheDeformedDomain,
                                                               testFactories,
                                                               testStiffnessVisitorsOnAnyDofsPerNode };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
