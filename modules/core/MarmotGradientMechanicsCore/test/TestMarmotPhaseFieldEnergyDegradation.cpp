#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot::Testing;
using namespace Marmot::NumericalAlgorithms::Differentiation;
using namespace Marmot::PhaseField::EnergyDegradationFunctions;

void testQuadratic()
{
  // Test quadratic degradation function g(pf) = (1-pf)^2
  // Boundary: g(0) = 1, g(1) = 0
  throwExceptionOnFailure( checkIfEqual( quadratic( 0.0 ), 1.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quadratic(0) != 1" );
  throwExceptionOnFailure( checkIfEqual( quadratic( 1.0 ), 0.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quadratic(1) != 0" );
  // Intermediate value: g(0.5) = 0.25
  throwExceptionOnFailure( checkIfEqual( quadratic( 0.5 ), 0.25 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quadratic(0.5) != 0.25" );
}

void testCubic()
{
  // Test cubic degradation function g(pf) = 3(1-pf)^2 - 2(1-pf)^3
  // Boundary: g(0) = 1, g(1) = 0
  throwExceptionOnFailure( checkIfEqual( cubic( 0.0 ), 1.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": cubic(0) != 1" );
  throwExceptionOnFailure( checkIfEqual( cubic( 1.0 ), 0.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": cubic(1) != 0" );
  // Intermediate value: g(0.5) = 0.5
  throwExceptionOnFailure( checkIfEqual( cubic( 0.5 ), 0.5 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": cubic(0.5) != 0.5" );
}

void testQuartic()
{
  // Test quartic degradation function g(pf) = 4(1-pf)^3 - 3(1-pf)^4
  // Boundary: g(0) = 1, g(1) = 0
  throwExceptionOnFailure( checkIfEqual( quartic( 0.0 ), 1.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quartic(0) != 1" );
  throwExceptionOnFailure( checkIfEqual( quartic( 1.0 ), 0.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quartic(1) != 0" );
  // Intermediate value: g(0.5) = 0.3125
  throwExceptionOnFailure( checkIfEqual( quartic( 0.5 ), 0.3125 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": quartic(0.5) != 0.3125" );
}

void testSecondOrderDerivedQuadratic()
{
  // SecondOrderDerived::quadratic returns {value, first_derivative, second_derivative}
  // g(pf) = (1-pf)^2,  dg/dpf = -2(1-pf),  d2g/dpf2 = 2

  // At pf = 0: g=1, dg=-2, d2g=2
  {
    auto [g, dg, d2g] = SecondOrderDerived::quadratic( 0.0 );
    throwExceptionOnFailure( checkIfEqual( g, 1.0 ), MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -2.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, 2.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0 failed" );
  }

  // At pf = 0.3: g=0.49, dg=-1.4, d2g=2
  {
    auto [g, dg, d2g] = SecondOrderDerived::quadratic( 0.3 );
    throwExceptionOnFailure( checkIfEqual( g, 0.49, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -1.4, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, 2.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0.3 failed" );
  }

  // At pf = 1: g=0, dg=0, d2g=2
  {
    auto [g, dg, d2g] = SecondOrderDerived::quadratic( 1.0 );
    throwExceptionOnFailure( checkIfEqual( g, 0.0 ), MakeString() << __PRETTY_FUNCTION__ << ": value at pf=1 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, 2.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=1 failed" );
  }

  // Consistency: returned value must match template function
  for ( const double pf : { 0.0, 0.25, 0.5, 0.75, 1.0 } ) {
    auto [g, dg, d2g] = SecondOrderDerived::quadratic( pf );
    throwExceptionOnFailure( checkIfEqual( g, quadratic( pf ) ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value inconsistency at pf=" << pf );
  }
}

void testSecondOrderDerivedCubic()
{
  // SecondOrderDerived::cubic returns {value, first_derivative, second_derivative}
  // g(pf) = 3s^2 - 2s^3  with s=1-pf
  // dg/dpf = -6s + 6s^2
  // d2g/dpf2 = 6 - 12s

  // At pf = 0.3: g=0.784, dg=-1.26, d2g=-2.4
  {
    auto [g, dg, d2g] = SecondOrderDerived::cubic( 0.3 );
    throwExceptionOnFailure( checkIfEqual( g, 0.784, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -1.26, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, -2.4, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0.3 failed" );
  }

  // At pf = 0.5: g=0.5, dg=-1.5, d2g=0
  {
    auto [g, dg, d2g] = SecondOrderDerived::cubic( 0.5 );
    throwExceptionOnFailure( checkIfEqual( g, 0.5 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0.5 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -1.5 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0.5 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, 0.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0.5 failed" );
  }

  // Consistency: returned value must match template function
  for ( const double pf : { 0.0, 0.25, 0.5, 0.75, 1.0 } ) {
    auto [g, dg, d2g] = SecondOrderDerived::cubic( pf );
    throwExceptionOnFailure( checkIfEqual( g, cubic( pf ) ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value inconsistency at pf=" << pf );
  }
}

void testSecondOrderDerivedQuartic()
{
  // SecondOrderDerived::quartic returns {value, first_derivative, second_derivative}
  // g(pf) = 4s^3 - 3s^4  with s=1-pf
  // dg/dpf = -12s^2 + 12s^3
  // d2g/dpf2 = 24s - 36s^2

  // At pf = 0.3: g=0.6517, dg=-1.764, d2g=-0.84
  {
    auto [g, dg, d2g] = SecondOrderDerived::quartic( 0.3 );
    throwExceptionOnFailure( checkIfEqual( g, 0.6517, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -1.764, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0.3 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, -0.84, 1e-14 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0.3 failed" );
  }

  // At pf = 0.5: g=0.3125, dg=-1.5, d2g=3
  {
    auto [g, dg, d2g] = SecondOrderDerived::quartic( 0.5 );
    throwExceptionOnFailure( checkIfEqual( g, 0.3125 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at pf=0.5 failed" );
    throwExceptionOnFailure( checkIfEqual( dg, -1.5 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative at pf=0.5 failed" );
    throwExceptionOnFailure( checkIfEqual( d2g, 3.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative at pf=0.5 failed" );
  }

  // Consistency: returned value must match template function
  for ( const double pf : { 0.0, 0.25, 0.5, 0.75, 1.0 } ) {
    auto [g, dg, d2g] = SecondOrderDerived::quartic( pf );
    throwExceptionOnFailure( checkIfEqual( g, quartic( pf ) ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value inconsistency at pf=" << pf );
  }
}

void testSecondOrderDerivedGeneric()
{
  // SecondOrderDerived::generic uses autodiff to compute derivatives
  // For p=2, a1=2, a2=0, a3=0: g(pf) = (1-pf)^2 / ((1-pf)^2 + 2*pf)
  // This is a known closed-form case

  // Boundary conditions must hold for all sensible parameter choices
  {
    auto [g0, dg0, d2g0] = SecondOrderDerived::generic( 0.0, 2.0, 2.0, 0.0, 0.0 );
    throwExceptionOnFailure( checkIfEqual( g0, 1.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": generic(pf=0) must equal 1" );
  }
  {
    auto [g1, dg1, d2g1] = SecondOrderDerived::generic( 1.0, 2.0, 2.0, 0.0, 0.0 );
    throwExceptionOnFailure( checkIfEqual( g1, 0.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": generic(pf=1) must equal 0" );
  }

  // Consistency: returned value must match template function
  for ( const double pf : { 0.0, 0.25, 0.5, 0.75 } ) {
    auto [g, dg, d2g] = SecondOrderDerived::generic( pf, 2.0, 2.0, 0.0, 0.0 );
    throwExceptionOnFailure( checkIfEqual( g, generic( pf, 2.0, 2.0, 0.0, 0.0 ) ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value inconsistency at pf=" << pf );
  }

  // Numerical derivative verification using central differences
  const double pf         = 0.4;
  const double p          = 2.0;
  const double a1         = 2.0;
  const double a2         = 0.5;
  const double a3         = 0.0;
  auto [g, dg, d2g]       = SecondOrderDerived::generic( pf, p, a1, a2, a3 );
  const double dg_numdiff = centralDifference( [p, a1, a2, a3](
                                                 const double x ) { return generic( x, p, a1, a2, a3 ); },
                                               pf );
  throwExceptionOnFailure( checkIfEqual( dg, dg_numdiff, 1e-6 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": first derivative mismatch (autodiff vs numdiff)" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testQuadratic,
                                                       testCubic,
                                                       testQuartic,
                                                       testSecondOrderDerivedQuadratic,
                                                       testSecondOrderDerivedCubic,
                                                       testSecondOrderDerivedQuartic,
                                                       testSecondOrderDerivedGeneric };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
