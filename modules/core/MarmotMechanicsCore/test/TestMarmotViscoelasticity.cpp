#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "autodiff/forward/real.hpp"

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Viscoelasticity;

void testAutodiffDerivative()
{
  // function parameters
  double m   = 0.5;
  double n   = 0.1;
  double tau = 1.234e5;

  // autodiff setup
  static constexpr int                         maxDerivativeOrder = 3;
  autodiff::Real< maxDerivativeOrder, double > tau_( tau );
  auto                                         phiPL = [&]( autodiff::Real< maxDerivativeOrder, double > t ) {
    return ComplianceFunctions::powerLaw( t, m, n );
  };
  auto phiLPL = [&]( autodiff::Real< maxDerivativeOrder, double > t ) {
    return ComplianceFunctions::logPowerLaw( t, m, n );
  };

  // autodiff derivatives
  Eigen::TensorFixedSize< double, Eigen::Sizes< maxDerivativeOrder > > autodiffValPL;
  Eigen::TensorFixedSize< double, Eigen::Sizes< maxDerivativeOrder > > autodiffValLPL;

  for ( int i = 0; i < maxDerivativeOrder; i++ ) {
    autodiffValPL( i )  = autodiff::derivatives( phiPL, autodiff::along( 1. ), autodiff::at( tau_ ) )[i + 1];
    autodiffValLPL( i ) = autodiff::derivatives( phiLPL, autodiff::along( 1. ), autodiff::at( tau_ ) )[i + 1];
  }

  // derivatives from sympy
  Eigen::TensorFixedSize< double, Eigen::Sizes< maxDerivativeOrder > > corrValPL;
  corrValPL.setValues( { m * n * std::pow( tau, n ) / tau,
                         m * n * std::pow( tau, n ) * ( n - 1 ) / std::pow( tau, 2 ),
                         m * n * std::pow( tau, n ) * ( std::pow( n, 2 ) - 3 * n + 2 ) / std::pow( tau, 3 ) } );

  Eigen::TensorFixedSize< double, Eigen::Sizes< maxDerivativeOrder > > corrValLPL;
  corrValLPL.setValues(
    { m * n * std::pow( tau, n ) / ( tau * ( std::pow( tau, n ) + 1 ) ),
      -m * n * std::pow( tau, n ) * ( n * std::pow( tau, n ) / ( std::pow( tau, n ) + 1 ) - n + 1 ) /
        ( std::pow( tau, 2 ) * ( std::pow( tau, n ) + 1 ) ),
      m * n * std::pow( tau, n ) *
        ( 2 * std::pow( n, 2 ) * std::pow( tau, 2 * n ) / std::pow( std::pow( tau, n ) + 1, 2 ) -
          3 * std::pow( n, 2 ) * std::pow( tau, n ) / ( std::pow( tau, n ) + 1 ) + std::pow( n, 2 ) +
          3 * n * std::pow( tau, n ) / ( std::pow( tau, n ) + 1 ) - 3 * n + 2 ) /
        ( std::pow( tau, 3 ) * ( std::pow( tau, n ) + 1 ) ) } );

  // check if equal
  throwExceptionOnFailure( checkIfEqual< double >( autodiffValPL, corrValPL ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " error in using autodiff on power law compliance function" );

  throwExceptionOnFailure( checkIfEqual< double >( autodiffValLPL, corrValLPL ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " error in using autodiff on log power law compliance function" );
}

int main()
{
  testAutodiffDerivative();
  return 0;
}
