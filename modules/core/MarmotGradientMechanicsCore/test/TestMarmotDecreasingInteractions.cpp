#include "Marmot/MarmotDecreasingInteractions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::NumericalAlgorithms::Differentiation;
using namespace Marmot::GradientDamage::DecreasingInteractions;

void testPohTemplate()
{
  // poh(omega, eta, R) = ((1-R)*exp(-eta*omega) + R - exp(-eta)) / (1 - exp(-eta))
  // At omega=0 (no damage): poh should equal 1 for any eta and R
  throwExceptionOnFailure( checkIfEqual( poh( 0.0, 1.0, 0.0 ), 1.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": poh(0, 1, 0) != 1" );
  throwExceptionOnFailure( checkIfEqual( poh( 0.0, 2.0, 0.3 ), 1.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": poh(0, 2, 0.3) != 1" );

  // At omega=1, eta=1, R=0: poh = (exp(-1) - exp(-1)) / (1 - exp(-1)) = 0
  throwExceptionOnFailure( checkIfEqual( poh( 1.0, 1.0, 0.0 ), 0.0 ),
                           MakeString() << __PRETTY_FUNCTION__ << ": poh(1, 1, 0) != 0" );

  // At omega=0.5, eta=1, R=0: expected value 0.37754066...
  const double expected = ( std::exp( -0.5 ) - std::exp( -1.0 ) ) / ( 1.0 - std::exp( -1.0 ) );
  throwExceptionOnFailure( checkIfEqual( poh( 0.5, 1.0, 0.0 ), expected ),
                           MakeString() << __PRETTY_FUNCTION__ << ": poh(0.5, 1, 0) mismatch" );

  // poh must be bounded in [R, 1]: poh(omega, eta, R) >= R for all omega in [0,1]
  for ( const double omega : { 0.1, 0.3, 0.5, 0.7, 0.9 } ) {
    const double val = poh( omega, 1.0, 0.1 );
    throwExceptionOnFailure( val >= 0.1 - 1e-14,
                             MakeString() << __PRETTY_FUNCTION__ << ": poh below R at omega=" << omega );
    throwExceptionOnFailure( val <= 1.0 + 1e-14,
                             MakeString() << __PRETTY_FUNCTION__ << ": poh above 1 at omega=" << omega );
  }
}

void testFirstOrderDerivedPoh()
{
  // FirstOrderDerived::poh returns {value, dvalue/domega}
  // For omega <= 0: returns {1, 0}
  {
    auto [g, dg] = FirstOrderDerived::poh( 0.0, 1.0, 0.0 );
    throwExceptionOnFailure( checkIfEqual( g, 1.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value at omega=0 must be 1" );
    throwExceptionOnFailure( checkIfEqual( dg, 0.0 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": derivative at omega=0 must be 0" );
  }

  // For omega > 0: value must match template function
  // dg/domega = -(1-R)*eta*exp(-eta*omega) / (1 - exp(-eta))
  {
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg]       = FirstOrderDerived::poh( omega, eta, R );

    const double expected_g  = poh( omega, eta, R );
    const double expected_dg = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                             MakeString() << __PRETTY_FUNCTION__ << ": derivative mismatch at omega=0.5" );
  }

  // Numerical derivative verification using central differences
  {
    const double omega      = 0.5;
    const double eta        = 1.0;
    const double R          = 0.0;
    auto [g, dg]            = FirstOrderDerived::poh( omega, eta, R );
    const double dg_numdiff = centralDifference( [eta, R]( const double w ) { return poh( w, eta, R ); }, omega );
    throwExceptionOnFailure( checkIfEqual( dg, dg_numdiff, 1e-6 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": derivative numdiff check failed at omega=0.5" );
  }
}

void testSecondOrderDerivedPoh()
{
  // SecondOrderDerived::poh returns {value, dvalue/domega, d2value/domega2}
  // d2g/domega2 = (1-R)*eta^2*exp(-eta*omega) / (1 - exp(-eta))
  // Note: -d2g/domega2 = -eta * dg/domega  (i.e. d2g = -eta * dg_domega)

  // At omega=0.5, eta=1, R=0: all three analytically known
  {
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg, d2g]  = SecondOrderDerived::poh( omega, eta, R );

    const double expected_g   = poh( omega, eta, R );
    const double expected_dg  = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    const double expected_d2g = ( 1.0 - R ) * eta * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );

    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                             MakeString() << __PRETTY_FUNCTION__ << ": first derivative mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( d2g, expected_d2g ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative mismatch at omega=0.5" );
  }

  // At omega=0.5, eta=2, R=0.2
  {
    const double omega = 0.5;
    const double eta   = 2.0;
    const double R     = 0.2;
    auto [g, dg, d2g]  = SecondOrderDerived::poh( omega, eta, R );

    const double expected_g   = poh( omega, eta, R );
    const double expected_dg  = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    const double expected_d2g = ( 1.0 - R ) * eta * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );

    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                             MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5, eta=2, R=0.2" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": first derivative mismatch at omega=0.5, eta=2, R=0.2" );
    throwExceptionOnFailure( checkIfEqual( d2g, expected_d2g ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << ": second derivative mismatch at omega=0.5, eta=2, R=0.2" );
  }

  // Numerical second derivative verification using central differences
  {
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg, d2g]  = SecondOrderDerived::poh( omega, eta, R );
    const double
      d2g_numdiff = centralDifference( [eta, R]( const double
                                                   w ) { return std::get< 1 >( FirstOrderDerived::poh( w, eta, R ) ); },
                                       omega );
    throwExceptionOnFailure( checkIfEqual( d2g, d2g_numdiff, 1e-6 ),
                             MakeString() << __PRETTY_FUNCTION__ << ": second derivative numdiff check failed" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testPohTemplate,
                                                       testFirstOrderDerivedPoh,
                                                       testSecondOrderDerivedPoh };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
