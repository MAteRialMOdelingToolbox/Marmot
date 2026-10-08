#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotWiechert.h"

#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::Materials::Wiechert;
using namespace Marmot::ContinuumMechanics::Viscoelasticity;

namespace {

  void testPowerLawRelaxationSpectrum()
  {
    constexpr int order   = 2;
    const double  m       = 12.0;
    const double  n       = 0.25;
    const double  minTau  = 0.01;
    const double  spacing = std::sqrt( 10. );

    const Properties relaxationTimes    = generateRelaxationTimes( 4, minTau, spacing );
    auto             relaxationFunction = [&]( autodiff::Real< order, double > time ) {
      return RelaxationFunctions::powerLaw( time, m, n );
    };

    const Properties elasticModuli = computeElasticModuli< order >( relaxationFunction, relaxationTimes, spacing );
    for ( int i = 0; i < relaxationTimes.size(); ++i ) {
      const double expectedSpectrum = m * n * ( n + 1. ) * std::pow( 2., -n ) * std::pow( relaxationTimes( i ), -n );
      throwExceptionOnFailure( checkIfEqual( evaluatePostWidderFormula< order >( relaxationFunction,
                                                                                 relaxationTimes( i ) ),
                                             expectedSpectrum,
                                             1e-12 ),
                               "Incorrect Post-Widder relaxation spectrum." );
      throwExceptionOnFailure( checkIfEqual( elasticModuli( i ), std::log( spacing ) * expectedSpectrum, 1e-12 ),
                               "Incorrect generalized-Maxwell branch modulus." );
      throwExceptionOnFailure( checkIfEqual( relaxationTimes( i ), minTau * std::pow( spacing, i ), 1e-14 ),
                               "Incorrect generalized-Maxwell relaxation time." );
    }
  }

  void testGaussQuadratureRelaxationSpectrum()
  {
    constexpr int order   = 2;
    const double  m       = 12.0;
    const double  n       = 0.25;
    const double  spacing = std::sqrt( 10. );

    const Properties relaxationTimes    = generateRelaxationTimes( 4, 0.01, spacing );
    auto             relaxationFunction = [&]( autodiff::Real< order, double > time ) {
      return RelaxationFunctions::powerLaw( time, m, n );
    };

    const Properties elasticModuli = computeElasticModuli< order >( relaxationFunction,
                                                                    relaxationTimes,
                                                                    spacing,
                                                                    /*gaussQuadrature*/ true );

    // two-point Gauss rule over one decade of the logarithmic relaxation-time axis
    const auto spectrum = [&]( double tau ) { return m * n * ( n + 1. ) * std::pow( 2., -n ) * std::pow( tau, -n ); };
    for ( int i = 0; i < relaxationTimes.size(); ++i ) {
      const double tau      = relaxationTimes( i );
      const double expected = std::log( spacing ) / 2. *
                              ( spectrum( tau * std::pow( spacing, -std::sqrt( 3. ) / 6. ) ) +
                                spectrum( tau * std::pow( spacing, std::sqrt( 3. ) / 6. ) ) );
      throwExceptionOnFailure( checkIfEqual( elasticModuli( i ), expected, 1e-12 ),
                               "Incorrect Gauss-quadrature generalized-Maxwell branch modulus." );
    }
  }

  void testStateVarMatrixIsUnchangedForVanishingTimeIncrement()
  {
    const Properties relaxationTimes = generateRelaxationTimes( 3, 0.01, std::sqrt( 10. ) );
    const Properties elasticModuli   = Properties::Constant( 3, 100.0 );

    StateVarMatrix   stateVars = StateVarMatrix::Ones( 6, 3 );
    Marmot::Vector6d dStrain;
    dStrain << 1e-3, 0., 0., 0., 0., 0.;

    updateStateVarMatrix( 0.0, elasticModuli, relaxationTimes, stateVars, dStrain, Marmot::Matrix6d::Identity() );
    throwExceptionOnFailure( stateVars.isOnes(), "State variables must not evolve when no time elapses." );

    updateStateVarMatrix( 1.0, elasticModuli, relaxationTimes, stateVars, dStrain, Marmot::Matrix6d::Identity() );
    throwExceptionOnFailure( !stateVars.isOnes(), "State variables must evolve over a finite time increment." );
  }

} // namespace

int main()
{
  executeTestsAndCollectExceptions( { testPowerLawRelaxationSpectrum,
                                      testGaussQuadratureRelaxationSpectrum,
                                      testStateVarMatrixIsUnchangedForVanishingTimeIncrement } );
  return 0;
}
