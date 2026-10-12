#include "Marmot/MarmotBulkViscosity.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <functional>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::FiniteElement::BulkViscosity;

namespace {
  // A representative concrete-like set, in the N/mm/tonne/s system the examples use.
  constexpr double density   = 2.4e-9;
  constexpr double waveSpeed = 3.7e6; // mm/s
  constexpr double length    = 12.5;  // mm
} // namespace

void testInactiveByDefault()
{
  const Coefficients defaultCoefficients;

  throwExceptionOnFailure( !defaultCoefficients.areActive(),
                           "Default-constructed bulk viscosity coefficients must be inactive." );

  // An inactive set must return EXACTLY zero, not merely something small: it is what keeps a run
  // that does not ask for bulk viscosity bit-identical to one built before it existed.
  for ( const double rate : { -1.0e3, -1.0, 0.0, 1.0, 1.0e3 } ) {
    const double stress = viscousStress( rate, density, waveSpeed, length, defaultCoefficients );
    throwExceptionOnFailure( stress == 0.0, "Inactive bulk viscosity must return exactly zero." );
  }
}

void testLinearTerm()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 0.0 };

  // Expansion: the linear term acts on both signs of the volumetric rate.
  const double rate     = 250.0;
  const double expected = 0.06 * density * waveSpeed * length * rate;

  throwExceptionOnFailure( checkIfEqual( viscousStress( rate, density, waveSpeed, length, coefficients ),
                                         expected,
                                         1e-14 ),
                           "Linear bulk viscosity term is wrong in expansion." );

  // Compression: same coefficient, opposite sign, and with the quadratic term switched off the
  // response must be exactly odd in the rate.
  throwExceptionOnFailure( checkIfEqual( viscousStress( -rate, density, waveSpeed, length, coefficients ),
                                         -expected,
                                         1e-14 ),
                           "Linear bulk viscosity term must be odd in the volumetric strain rate." );
}

void testQuadraticTermIsCompressionOnly()
{
  const Coefficients linearOnly{ .linear = 0.06, .quadratic = 0.0 };
  const Coefficients both{ .linear = 0.06, .quadratic = 1.2 };

  const double rate = 250.0;

  // In expansion the quadratic term must contribute nothing at all, so that an opening crack is
  // not resisted by an artificial stress.
  throwExceptionOnFailure( viscousStress( rate, density, waveSpeed, length, both ) ==
                             viscousStress( rate, density, waveSpeed, length, linearOnly ),
                           "The quadratic term must be inactive in expansion." );

  // In compression it adds a compressive (negative) contribution of the documented magnitude.
  const double lengthTimesB2 = 1.2 * length;
  const double expected      = 0.06 * density * waveSpeed * length * ( -rate ) -
                          density * lengthTimesB2 * lengthTimesB2 * rate * rate;

  throwExceptionOnFailure( checkIfEqual( viscousStress( -rate, density, waveSpeed, length, both ), expected, 1e-14 ),
                           "Quadratic bulk viscosity term is wrong in compression." );
}

void testDissipativity()
{
  // The whole point of the term is that it can only ever REMOVE energy. The power density is
  // sigma_bv * rate, and it must be non-negative for every rate and every admissible coefficient
  // pair -- otherwise the device that is meant to damp the solution could drive it instead.
  const std::vector< Coefficients > coefficientSets = { { .linear = 0.06, .quadratic = 0.0 },
                                                        { .linear = 0.0, .quadratic = 1.2 },
                                                        { .linear = 0.06, .quadratic = 1.2 },
                                                        { .linear = 2.0, .quadratic = 5.0 } };

  for ( const auto& coefficients : coefficientSets ) {
    for ( double rate = -1.0e4; rate <= 1.0e4; rate += 137.0 ) {
      const double power = viscousStress( rate, density, waveSpeed, length, coefficients ) * rate;
      throwExceptionOnFailure( power >= 0.0, "Bulk viscosity must never do positive work on the system." );
    }
  }
}

void testZeroRateGivesZeroStress()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 1.2 };

  throwExceptionOnFailure( viscousStress( 0.0, density, waveSpeed, length, coefficients ) == 0.0,
                           "A vanishing volumetric strain rate must produce no viscous stress." );
}

void testScaling()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 0.0 };
  const double       rate = -400.0;

  const double reference = viscousStress( rate, density, waveSpeed, length, coefficients );

  // The linear term is first order in each of rho, c and L.
  throwExceptionOnFailure( checkIfEqual( viscousStress( rate, 2.0 * density, waveSpeed, length, coefficients ),
                                         2.0 * reference,
                                         1e-14 ),
                           "The linear term must scale linearly with the density." );
  throwExceptionOnFailure( checkIfEqual( viscousStress( rate, density, 2.0 * waveSpeed, length, coefficients ),
                                         2.0 * reference,
                                         1e-14 ),
                           "The linear term must scale linearly with the wave speed." );
  throwExceptionOnFailure( checkIfEqual( viscousStress( rate, density, waveSpeed, 2.0 * length, coefficients ),
                                         2.0 * reference,
                                         1e-14 ),
                           "The linear term must scale linearly with the characteristic length." );

  // The quadratic term is second order in L.
  const Coefficients quadraticOnly{ .linear = 0.0, .quadratic = 1.2 };
  const double       quadraticReference = viscousStress( rate, density, waveSpeed, length, quadraticOnly );
  throwExceptionOnFailure( checkIfEqual( viscousStress( rate, density, waveSpeed, 2.0 * length, quadraticOnly ),
                                         4.0 * quadraticReference,
                                         1e-14 ),
                           "The quadratic term must scale with the square of the characteristic length." );
}

void testIncrementOverload()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 1.2 };

  const double timeIncrement = 5.0e-8;
  const double increment     = -1.25e-5;
  const double rate          = increment / timeIncrement;

  throwExceptionOnFailure( checkIfEqual( viscousStressFromIncrement( increment,
                                                                     timeIncrement,
                                                                     density,
                                                                     waveSpeed,
                                                                     length,
                                                                     coefficients ),
                                         viscousStress( rate, density, waveSpeed, length, coefficients ),
                                         1e-14 ),
                           "The increment overload must agree with the rate form." );

  // An explicit solver primes its internal force with dT = 0. That must not divide by zero.
  for ( const double nonPositive : { 0.0, -1.0e-8 } ) {
    const double stress = viscousStressFromIncrement( increment,
                                                      nonPositive,
                                                      density,
                                                      waveSpeed,
                                                      length,
                                                      coefficients );
    throwExceptionOnFailure( stress == 0.0 && std::isfinite( stress ),
                             "A non-positive time increment must yield exactly zero viscous stress." );
  }
}

void testDegradationIsOffByDefault()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 1.2 };

  throwExceptionOnFailure( !coefficients.isDegraded(),
                           "The degradation with damage must be off unless it is asked for." );

  // A zero exponent must not even look at the wave speeds: it is the branch that keeps a deck
  // which does not request the degradation bit-identical to one predating it.
  throwExceptionOnFailure( degradationFactor( 0.0, waveSpeed, coefficients.degradation ) == 1.0,
                           "A zero exponent must return exactly one, whatever the wave speeds are." );
}

void testDegradationExponents()
{
  // Half the wave speed is a quarter of the tangent stiffness, so for a model whose tangent
  // degrades as (1 - omega) this is omega = 0.75.
  const double current = 0.5 * waveSpeed;

  throwExceptionOnFailure( checkIfEqual( degradationFactor( current, waveSpeed, 1.0 ), 0.5, 1e-14 ),
                           "Exponent one must scale with the wave speed itself, i.e. sqrt(1 - omega)." );
  throwExceptionOnFailure( checkIfEqual( degradationFactor( current, waveSpeed, 2.0 ), 0.25, 1e-14 ),
                           "Exponent two must scale with the tangent stiffness, i.e. (1 - omega)." );

  // A non-integral exponent goes through std::pow; check it against the closed form so that the
  // fast paths for 1 and 2 cannot drift away from the general one.
  throwExceptionOnFailure( checkIfEqual( degradationFactor( current, waveSpeed, 0.5 ), std::sqrt( 0.5 ), 1e-14 ),
                           "A non-integral exponent must agree with the closed form." );
}

void testDegradationIsClampedAndSafe()
{
  throwExceptionOnFailure( degradationFactor( 0.0, waveSpeed, 2.0 ) == 0.0,
                           "A fully degraded material must carry no artificial viscous stress: that is the whole "
                           "point, since the linear term would otherwise resist the crack opening." );

  // A tangent stiffer than the reference must not AMPLIFY the damping above the value the
  // undamaged reference sets -- compaction and hardening can both do that.
  throwExceptionOnFailure( degradationFactor( 3.0 * waveSpeed, waveSpeed, 1.0 ) == 1.0,
                           "The factor must be clamped at one, never amplifying the viscous stress." );

  // A reference of zero is a material that was never asked, or a massless one; dividing by it
  // would be a silent infinity.
  throwExceptionOnFailure( degradationFactor( waveSpeed, 0.0, 2.0 ) == 1.0,
                           "A non-positive reference wave speed must fall back to no degradation." );
}

void testDegradationScalesTheViscousStress()
{
  const Coefficients coefficients{ .linear = 0.06, .quadratic = 1.2, .degradation = 2.0 };
  const double       rate = 400.0;

  throwExceptionOnFailure( coefficients.isDegraded(), "A positive exponent must report the term as degraded." );

  // What the element multiplies: the degradation is a factor on the whole viscous stress, so an
  // element at omega = 0.75 integrates a quarter of it.
  const double undegraded = viscousStress( rate, density, waveSpeed, length, coefficients );
  const double factor     = degradationFactor( 0.5 * waveSpeed, waveSpeed, coefficients.degradation );

  throwExceptionOnFailure( checkIfEqual( factor * undegraded, 0.25 * undegraded, 1e-14 ),
                           "The degraded viscous stress must be the factor times the undegraded one." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testInactiveByDefault,
                                                       testLinearTerm,
                                                       testQuadraticTermIsCompressionOnly,
                                                       testDissipativity,
                                                       testZeroRateGivesZeroStress,
                                                       testScaling,
                                                       testIncrementOverload,
                                                       testDegradationIsOffByDefault,
                                                       testDegradationExponents,
                                                       testDegradationIsClampedAndSafe,
                                                       testDegradationScalesTheViscousStress };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
