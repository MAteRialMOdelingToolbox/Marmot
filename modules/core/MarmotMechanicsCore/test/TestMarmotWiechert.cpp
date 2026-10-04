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

} // namespace

int main()
{
  executeTestsAndCollectExceptions( { testPowerLawRelaxationSpectrum } );
  return 0;
}
