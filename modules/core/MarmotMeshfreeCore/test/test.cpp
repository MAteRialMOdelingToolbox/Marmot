#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed.h"
#include "Marmot/MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximation.h"
#include "Marmot/MarmotMeshfreeReproducingKernelApproximationImplicit.h"
#include "Marmot/MarmotMonomialBasisFunctions.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Meshfree;
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
    std::vector< Eigen::VectorXd >                     centers;
    std::vector< std::unique_ptr< Kernel > >           kernels;
    std::vector< const MarmotMeshfreeKernelFunction* > pointers;

    KernelGrid( int dim, int n, double supportRadius )
    {
      const int nTotal = std::pow( n, dim );
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

} // namespace

void testMonomialBasis()
{
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 0; order <= 3; order++ ) {
      const auto a = exponents( order, dim );
      const int  n = Math::computeSizeOfMonomialBasisVector( order, dim );
      throwExceptionOnFailure( n == int( a.size() ),
                               MakeString() << "size of the basis of order " << order << " in " << dim << "D" );

      const Eigen::VectorXd x = point( dim, 0.7 );
      Eigen::VectorXd       H( n );
      Eigen::MatrixXd       dH( n, dim );
      Math::computeMonomialBasis( order, x, H );
      Math::computeMonomialBasisGradient( order, x, dH );

      // the basis is ordered as the exponent tuples above; compare with the direct evaluation
      for ( int k = 0; k < n; k++ ) {
        throwExceptionOnFailure( checkIfEqual( H[k], monomial( a[k], x ), 1e-14 ),
                                 MakeString() << "monomial " << k << " of order " << order << " in " << dim << "D" );
        for ( int i = 0; i < dim; i++ )
          throwExceptionOnFailure( checkIfEqual( dH( k, i ), monomialDerivative( a[k], x, i ), 1e-13 ),
                                   MakeString()
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
                                 MakeString() << name << ": gradient " << i << " in " << dim << "D" );
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
  checkKernel< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed >( "BSpline2ndOrderBoxed" );
  checkKernel< MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed >( "BSpline3rdOrderBoxed" );
}

// reproduction of the complete polynomial basis of the requested order by the values, and of its derivatives by
// the gradients; for the explicit (direct) derivatives, also consistency with the shape function values
template < typename Kernel, typename Approximation >
void checkApproximation( const std::string& name, int dim, int order, double supportRadius, bool directDerivatives )
{
  KernelGrid< Kernel > grid( dim, 6, supportRadius );
  const Approximation  approximation( dim, order );
  const int            nNodes = grid.pointers.size();
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
                               MakeString()
                                 << name << ", order " << order << " in " << dim << "D: monomial not reproduced" );
      for ( int i = 0; i < dim; i++ )
        throwExceptionOnFailure( std::abs( reproducedGradient[i] - monomialDerivative( alpha, x, i ) ) < 1e-9,
                                 MakeString() << name << ", order " << order << " in " << dim
                                              << "D: derivative of a monomial not reproduced, " << reproducedGradient[i]
                                              << " vs " << monomialDerivative( alpha, x, i ) );
    }

    if ( directDerivatives )
      for ( int i = 0; i < dim; i++ ) {
        const double    h  = 1e-6;
        Eigen::VectorXd xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        const Eigen::VectorXd fd = ( values( xp ) - values( xm ) ) / ( 2 * h );
        throwExceptionOnFailure( ( dN.row( i ).transpose() - fd ).norm() < 1e-6 * ( 1 + fd.norm() ),
                                 MakeString() << name << ", order " << order << " in " << dim
                                              << "D: gradients are not the derivatives of the values" );
      }
  }

  Eigen::MatrixXd dN( dim, nNodes );
  bool            thrown = false;
  try {
    approximation.computeShapeFunctionGradients( point( dim, 2.1 ).data(), grid.pointers, dN.data() );
  }
  catch ( const std::runtime_error& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, name + ": computeShapeFunctionGradients is not implemented and must throw" );
}

void testReproducingKernelApproximation()
{
  using K2 = MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed;
  using K3 = MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed;
  using RK = MarmotMeshfreeReproducingKernelApproximation;
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      const double r = order == 1 ? 1.6 : 2.4;
      checkApproximation< K2, RK >( "RK, 2nd order B-spline", dim, order, r, true );
      checkApproximation< K3, RK >( "RK, 3rd order B-spline", dim, order, r, true );
    }
}

void testImplicitGradientReproducingKernelApproximation()
{
  using K2  = MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed;
  using K3  = MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed;
  using IRK = MarmotMeshfreeReproducingKernelApproximationImplicit;
  // implicit gradients reproduce the derivatives of the basis, but are not the derivatives of the values
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      const double r = order == 1 ? 1.6 : 2.4;
      checkApproximation< K2, IRK >( "implicit RK, 2nd order B-spline", dim, order, r, false );
      checkApproximation< K3, IRK >( "implicit RK, 3rd order B-spline", dim, order, r, false );
    }
}

// the moment matrix helpers are protected: expose them for the test
struct MomentMatrixAccess : MarmotMeshfreeReproducingKernelApproximation {
  using MarmotMeshfreeReproducingKernelApproximation::computeMMatrix;
  using MarmotMeshfreeReproducingKernelApproximation::computeMMatrixAndGradient;
};

void testMomentMatrixGradient()
{
  using RK = MomentMatrixAccess;
  for ( int dim = 1; dim <= 3; dim++ )
    for ( int order = 1; order <= 2; order++ ) {
      KernelGrid< MarmotMeshfreeKernelFunctionBSpline3rdOrderBoxed > grid( dim, 6, 2.4 );
      const Eigen::VectorXd                                          x = point( dim, 2.3 );
      const auto [M, dM] = RK::computeMMatrixAndGradient( x, grid.pointers, order );
      throwExceptionOnFailure( ( M - RK::computeMMatrix( x, grid.pointers, order ) ).norm() < 1e-14 * M.norm(),
                               "moment matrix of both evaluations differs" );
      for ( int i = 0; i < dim; i++ ) {
        const double    h  = 1e-6;
        Eigen::VectorXd xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        const Eigen::MatrixXd fd = ( RK::computeMMatrix( xp, grid.pointers, order ) -
                                     RK::computeMMatrix( xm, grid.pointers, order ) ) /
                                   ( 2 * h );
        throwExceptionOnFailure( ( dM[i] - fd ).norm() < 1e-6 * ( 1 + fd.norm() ),
                                 MakeString() << "gradient " << i << " of the moment matrix of order " << order
                                              << " in " << dim << "D" );
      }
    }
}

void testCompletenessOrderIsReducedForFewNodes()
{
  // two kernels in 1D cannot carry a quadratic basis: the order drops to 1
  KernelGrid< MarmotMeshfreeKernelFunctionBSpline2ndOrderBoxed > grid( 1, 2, 1.6 );
  const MarmotMeshfreeReproducingKernelApproximation             approximation( 1, 2 );
  const double                                                   x = 0.4;
  double                                                         N[2], dN[2];
  approximation.computeShapeFunctionsAndGradients( &x, grid.pointers, N, dN );
  throwExceptionOnFailure( checkIfEqual( N[0] + N[1], 1.0, 1e-14 ), "partition of unity with a reduced order" );
  throwExceptionOnFailure( std::abs( N[0] * 0 + N[1] * 1 - x ) < 1e-12,
                           "two nodes in 1D still reproduce a linear field" );
}

int main()
{
  auto testFunctions = std::vector< std::function< void() > >{ testMonomialBasis,
                                                               testKernelFunctions,
                                                               testReproducingKernelApproximation,
                                                               testImplicitGradientReproducingKernelApproximation,
                                                               testMomentMatrixGradient,
                                                               testCompletenessOrderIsReducedForFewNodes };
  executeTestsAndCollectExceptions( testFunctions );
  return 0;
}
