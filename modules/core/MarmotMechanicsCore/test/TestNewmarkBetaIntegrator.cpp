#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <cmath>

using namespace Marmot::Testing;

using namespace Marmot;

using Marmot::ContinuumMechanics::TimeIntegration::newmarkBetaIntegration;

/**
 * @brief Checks a value against its expectation, naming the calling test and the quantity in the failure message.
 * @param[in] what The quantity.
 * @param[in] value Its value.
 * @param[in] expected The expected value.
 * @param[in] caller The calling test.
 */
void checkValue( const char* what, double value, double expected, const char* caller )
{
  throwExceptionOnFailure( checkIfEqual( value, expected, 1e-10 ),
                           Marmot::Utilities::MakeString()
                             << caller << " failed: " << what << " = " << value << ", expected " << expected );
}

/**
 * @brief The average acceleration rule from rest: a unit increment over a unit time step gives a = 2 du / dT^2 and
 * v = a dT / 2.
 */
void test_NewmarkBetaIntegrator()
{
  const double beta  = 0.25;
  const double gamma = 0.5;

  double du = 0.5;
  double dT = 1.0;

  double v     = 0;
  double a     = 0;
  double da_du = 0;

  newmarkBetaIntegration< 1 >( &du, &v, &a, dT, beta, gamma, &da_du );

  checkValue( "a", a, 2, __PRETTY_FUNCTION__ );
  checkValue( "v", v, 1, __PRETTY_FUNCTION__ );
  checkValue( "da/ddu", da_du, 1. / ( beta * dT * dT ), __PRETTY_FUNCTION__ );
}

/**
 * @brief A zero-length increment adds no inertia: velocity and acceleration are left as they are, whatever the
 * displacement increment, and the derivative of the acceleration is zero. The former clamp of the time increment to
 * 1e-16 instead produced an acceleration of order -v / (beta 1e-16).
 */
void test_NewmarkBetaIntegratorZeroIncrement()
{
  const double beta  = 0.25;
  const double gamma = 0.5;

  for ( const double du : { 0.0, 0.3 } ) {
    double       v[2]     = { 1.5, -2.0 };
    double       a[2]     = { 0.25, 4.0 };
    double       da_du[4] = { -1.0, -1.0, -1.0, -1.0 };
    double       duVec[2] = { du, -du };
    const double dT0      = 0.0;

    newmarkBetaIntegration< 2 >( duVec, v, a, dT0, beta, gamma, da_du );

    checkValue( "v[0]", v[0], 1.5, __PRETTY_FUNCTION__ );
    checkValue( "v[1]", v[1], -2.0, __PRETTY_FUNCTION__ );
    checkValue( "a[0]", a[0], 0.25, __PRETTY_FUNCTION__ );
    checkValue( "a[1]", a[1], 4.0, __PRETTY_FUNCTION__ );
    checkValue( "da/ddu[0,0]", da_du[0], 0.0, __PRETTY_FUNCTION__ );
    checkValue( "da/ddu[1,1]", da_du[3], 0.0, __PRETTY_FUNCTION__ );
    // the off-diagonal entries are not written
    checkValue( "da/ddu[1,0]", da_du[1], -1.0, __PRETTY_FUNCTION__ );
  }
}

/**
 * @brief beta = 0 switches the time integration off: zero acceleration and derivative, and the velocity
 * v_n + dT (1 - gamma) a_n.
 */
void test_NewmarkBetaIntegratorQuasiStatic()
{
  double du    = 0.5;
  double v     = 1.0;
  double a     = 2.0;
  double da_du = -1.0;

  newmarkBetaIntegration< 1 >( &du, &v, &a, 0.5, 0.0, 0.0, &da_du );

  checkValue( "a", a, 0.0, __PRETTY_FUNCTION__ );
  checkValue( "v", v, 2.0, __PRETTY_FUNCTION__ );
  checkValue( "da/ddu", da_du, 0.0, __PRETTY_FUNCTION__ );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ test_NewmarkBetaIntegrator,
                                                       test_NewmarkBetaIntegratorZeroIncrement,
                                                       test_NewmarkBetaIntegratorQuasiStatic };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
