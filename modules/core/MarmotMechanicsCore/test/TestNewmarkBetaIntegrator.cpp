#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/NewmarkBetaIntegrator.h"
#include <cmath>

using namespace Marmot::Testing;

using namespace Marmot;

void test_NewmarkBetaIntegrator()
{
  const double beta  = 0.25;
  const double gamma = 0.5;

  double du = 0.5;
  double dT = 1.0;

  double v     = 0;
  double a     = 0;
  double da_du = 0;

  TimeIntegration::newmarkBetaIntegration< 1 >( &du, &v, &a, dT, beta, gamma, &da_du );

  throwExceptionOnFailure( checkIfEqual( a, 2, 1e-10 ), MakeString() << __PRETTY_FUNCTION__ << " failed: " << a );
  throwExceptionOnFailure( checkIfEqual( v, 1, 1e-10 ), MakeString() << __PRETTY_FUNCTION__ << " failed: " << v );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ test_NewmarkBetaIntegrator };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
