#include "Marmot/MarmotIsotropicTensorFunctions.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Math;
using namespace Marmot::FastorStandardTensors;

// ---------------------------------------------------------------------------
// Test matrices: distinct eigenvalues, fully repeated (A=cI, the critical case --
// every reference/undeformed configuration in a hyperelastic potential), and
// partially repeated (two equal, one distinct).
// ---------------------------------------------------------------------------

Tensor33d distinctA()
{
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 2.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 5.0;
  A( 0, 1 ) = A( 1, 0 ) = 0.3;
  A( 0, 2 ) = A( 2, 0 ) = -0.2;
  A( 1, 2 ) = A( 2, 1 ) = 0.1;
  return A;
}

Tensor33d fullyRepeatedA()
{
  Tensor33d A( 0.0 );
  A( 0, 0 ) = A( 1, 1 ) = A( 2, 2 ) = 4.0;
  return A;
}

Tensor33d partiallyRepeatedA()
{
  Tensor33d A( 0.0 );
  A( 0, 0 ) = 3.0;
  A( 1, 1 ) = 3.0;
  A( 2, 2 ) = 7.0;
  return A;
}

std::vector< std::pair< std::string, Tensor33d > > testMatrices()
{
  return { { "distinct", distinctA() },
           { "fully repeated (A=cI)", fullyRepeatedA() },
           { "partially repeated", partiallyRepeatedA() } };
}

Tensor33d symmetricallyPerturbed( Tensor33d A, int i, int j, double delta )
{
  A( i, j ) += delta;
  if ( i != j )
    A( j, i ) += delta;
  return A;
}

// ---------------------------------------------------------------------------
// P-1/P-2: first derivative matches known closed forms h(x)=x (g=tr(A), dg=I,
// d2g=0) and h(x)=x^2 (g=tr(A^2), dg=2A, d2g=2*I4sym), at every eigenvalue
// multiplicity -- including at A=cI, the exact configuration where a generic
// AD-through-eigendecomposition approach silently discards derivative content.
// ---------------------------------------------------------------------------

void testKnownClosedForms()
{
  auto h1   = []( double x ) { return x; };
  auto h1p  = []( double ) { return 1.0; };
  auto h1pp = []( double ) { return 0.0; };

  auto h2   = []( double x ) { return x * x; };
  auto h2p  = []( double x ) { return 2.0 * x; };
  auto h2pp = []( double ) { return 2.0; };

  Tensor33d   I3 = Spatial3D::I;
  Tensor3333d zero4( 0.0 );
  using namespace Fastor;
  using namespace Marmot::FastorIndices;
  Tensor3333d I4sym = 0.5 * ( einsum< ik, jl, to_ijkl >( I3, I3 ) + einsum< il, jk, to_ijkl >( I3, I3 ) );

  for ( const auto& [name, A] : testMatrices() ) {
    {
      auto [g, dg, d2g] = sumOfEigenvaluesFunctionAndDerivatives( A, h1, h1p, h1pp );
      throwExceptionOnFailure( checkIfEqual( g, Fastor::trace( A ), 1e-10 ),
                               "P-1 [" + name + "]: h=x energy does not match trace(A)" );
      throwExceptionOnFailure( checkIfEqual( dg, I3, 1e-10 ),
                               "P-1 [" + name + "]: h=x first derivative does not match I" );
      throwExceptionOnFailure( checkIfEqual( d2g, zero4, 1e-10 ),
                               "P-1 [" + name + "]: h=x second derivative does not match 0" );
    }
    {
      auto [g, dg, d2g]    = sumOfEigenvaluesFunctionAndDerivatives( A, h2, h2p, h2pp );
      const double trA2    = Fastor::trace( Fastor::matmul( A, A ) );
      Tensor33d    dgKnown = 2.0 * A;
      throwExceptionOnFailure( checkIfEqual( g, trA2, 1e-10 ),
                               "P-2 [" + name + "]: h=x^2 energy does not match trace(A^2)" );
      throwExceptionOnFailure( checkIfEqual( dg, dgKnown, 1e-10 ),
                               "P-2 [" + name + "]: h=x^2 first derivative does not match 2*A" );
      Tensor3333d d2gKnown = 2.0 * I4sym;
      throwExceptionOnFailure( checkIfEqual( d2g, d2gKnown, 1e-10 ),
                               "P-2 [" + name + "]: h=x^2 second derivative does not match 2*I4sym" );
    }
  }
}

// ---------------------------------------------------------------------------
// I-1: first derivative matches a raw (independent) finite difference of the
// scalar g itself, for a non-integer exponent mimicking a real Ogden term.
// ---------------------------------------------------------------------------

void testFirstDerivativeMatchesFiniteDifference()
{
  auto h   = []( double x ) { return std::pow( x, 2.5 ); };
  auto hp  = []( double x ) { return 2.5 * std::pow( x, 1.5 ); };
  auto hpp = []( double x ) { return 2.5 * 1.5 * std::pow( x, 0.5 ); };

  const double fdH = 1e-6;

  for ( const auto& [name, A] : testMatrices() ) {
    auto [g, dg, d2g] = sumOfEigenvaluesFunctionAndDerivatives( A, h, hp, hpp );

    Tensor33d fdDg( 0.0 );
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        const double delta = ( i == j ) ? fdH : fdH / 2.0;
        auto [gp,
              _1,
              _2]    = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( A, i, j, delta ), h, hp, hpp );
        auto [gm,
              _3,
              _4]    = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( A, i, j, -delta ), h, hp, hpp );
        fdDg( i, j ) = ( gp - gm ) / ( 2.0 * fdH );
      }
    }

    throwExceptionOnFailure( checkIfEqual( dg, fdDg, 1e-4 ),
                             "I-1 [" + name + "]: h=x^2.5 first derivative does not match raw finite difference" );
  }
}

// ---------------------------------------------------------------------------
// I-2: second derivative matches a raw finite difference of the (already
// validated) FIRST derivative -- an independent check of the second-derivative
// closed form (not comparing the formula to itself), for the same non-integer
// exponent, at every eigenvalue multiplicity.
// ---------------------------------------------------------------------------

void testSecondDerivativeMatchesFiniteDifferenceOfFirst()
{
  auto h   = []( double x ) { return std::pow( x, 2.5 ); };
  auto hp  = []( double x ) { return 2.5 * std::pow( x, 1.5 ); };
  auto hpp = []( double x ) { return 2.5 * 1.5 * std::pow( x, 0.5 ); };

  const double fdH = 1e-6;

  for ( const auto& [name, A] : testMatrices() ) {
    auto [g, dg, d2g] = sumOfEigenvaluesFunctionAndDerivatives( A, h, hp, hpp );

    Tensor3333d fdD2g( 0.0 );
    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        const double delta = ( k == l ) ? fdH : fdH / 2.0;
        auto [_1,
              dgp,
              _2] = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( A, k, l, delta ), h, hp, hpp );
        auto [_3,
              dgm,
              _4] = sumOfEigenvaluesFunctionAndDerivatives( symmetricallyPerturbed( A, k, l, -delta ), h, hp, hpp );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            fdD2g( i, j, k, l ) = ( dgp( i, j ) - dgm( i, j ) ) / ( 2.0 * fdH );
      }
    }

    throwExceptionOnFailure( checkIfEqual( d2g, fdD2g, 1e-3 ),
                             "I-2 [" + name +
                               "]: h=x^2.5 second derivative does not match finite difference of first derivative" );
  }
}

// ---------------------------------------------------------------------------
// I-3: third derivative (from sumOfEigenvaluesFunctionAndDerivativesUpToThird,
// itself computed as a finite difference of the closed-form second derivative)
// matches an INDEPENDENT raw directional triple finite difference of g itself
// -- the standard 8-corner stencil, bypassing this module's own first/second
// derivative machinery entirely. Direction triples deliberately include PURE
// OFF-DIAGONAL perturbations at A=cI: exactly the configuration that silently
// broke under differentiating computeEigenSystemJacobi directly.
// ---------------------------------------------------------------------------

double rawDirectionalThirdFiniteDifference( const Tensor33d&                         A,
                                            const Tensor33d&                         dA1,
                                            const Tensor33d&                         dA2,
                                            const Tensor33d&                         dA3,
                                            const std::function< double( double ) >& h,
                                            const std::function< double( double ) >& hp,
                                            const std::function< double( double ) >& hpp,
                                            double                                   s = 1e-3 )
{
  const double sign[2] = { 1.0, -1.0 };
  double       total   = 0.0;
  for ( int a = 0; a < 2; ++a ) {
    for ( int b = 0; b < 2; ++b ) {
      for ( int c = 0; c < 2; ++c ) {
        Tensor33d Apert     = A + sign[a] * s * dA1 + sign[b] * s * dA2 + sign[c] * s * dA3;
        auto [gval, _1, _2] = sumOfEigenvaluesFunctionAndDerivatives( Apert, h, hp, hpp );
        total += sign[a] * sign[b] * sign[c] * gval;
      }
    }
  }
  return total / ( 8.0 * s * s * s );
}

double contractThirdDerivative( const Tensor333333d& d3g,
                                const Tensor33d&     dA1,
                                const Tensor33d&     dA2,
                                const Tensor33d&     dA3 )
{
  double total = 0.0;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          for ( int m = 0; m < 3; ++m )
            for ( int n = 0; n < 3; ++n )
              total += d3g( i, j, k, l, m, n ) * dA1( i, j ) * dA2( k, l ) * dA3( m, n );
  return total;
}

void testThirdDerivativeMatchesIndependentRawFiniteDifference()
{
  auto h   = []( double x ) { return std::pow( x, 2.5 ); };
  auto hp  = []( double x ) { return 2.5 * std::pow( x, 1.5 ); };
  auto hpp = []( double x ) { return 2.5 * 1.5 * std::pow( x, 0.5 ); };

  Tensor33d E01( 0.0 );
  E01( 0, 1 ) = E01( 1, 0 ) = 1.0;
  Tensor33d E02( 0.0 );
  E02( 0, 2 ) = E02( 2, 0 ) = 1.0;
  Tensor33d E12( 0.0 );
  E12( 1, 2 ) = E12( 2, 1 ) = 1.0;
  Tensor33d E00( 0.0 );
  E00( 0, 0 ) = 1.0;
  Tensor33d E11( 0.0 );
  E11( 1, 1 ) = 1.0;
  Tensor33d Emix( 0.0 );
  Emix( 0, 0 ) = 0.3;
  Emix( 1, 1 ) = -0.4;
  Emix( 0, 1 ) = Emix( 1, 0 ) = 0.5;
  Emix( 0, 2 ) = Emix( 2, 0 ) = -0.2;

  struct Dirs {
    Tensor33d   d1, d2, d3;
    std::string label;
  };
  std::vector< Dirs > dirSets = { { E01, E02, E12, "pure off-diagonal triple" },
                                  { E00, E11, E01, "diag+diag+offdiag" },
                                  { Emix, E01, E02, "generic mixed" } };

  for ( const auto& [name, A] : testMatrices() ) {
    auto [g, dg, d2g, d3g] = sumOfEigenvaluesFunctionAndDerivativesUpToThird( A, h, hp, hpp );

    for ( const auto& d : dirSets ) {
      const double analytic = contractThirdDerivative( d3g, d.d1, d.d2, d.d3 );
      const double raw      = rawDirectionalThirdFiniteDifference( A, d.d1, d.d2, d.d3, h, hp, hpp );
      const double relErr   = std::abs( analytic - raw ) / ( 1.0 + std::abs( raw ) );
      throwExceptionOnFailure( relErr < 1e-3,
                               "I-3 [" + name + " / " + d.label +
                                 "]: third derivative does not match independent raw directional finite "
                                 "difference (analytic=" +
                                 std::to_string( analytic ) + ", raw=" + std::to_string( raw ) + ")" );
    }
  }
}

int main()
{
  std::vector< std::function< void() > > tests = {
    testKnownClosedForms,
    testFirstDerivativeMatchesFiniteDifference,
    testSecondDerivativeMatchesFiniteDifferenceOfFirst,
    testThirdDerivativeMatchesIndependentRawFiniteDifference,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
