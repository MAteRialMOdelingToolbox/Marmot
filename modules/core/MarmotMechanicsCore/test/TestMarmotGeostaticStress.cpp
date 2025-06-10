#include "Marmot/MarmotGeostaticStress.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <cmath>

using namespace Marmot::Testing;

using namespace Marmot;

void test_geostaticStress()
{

  const double geostaticStressDefintion[] = { 0, 0, 200.0, 100.0, 0.5, 0.2 };

  const auto [s11, s22, s33] = GeostaticStress::getGeostaticStressFromLinearDistribution( geostaticStressDefintion,
                                                                                          50 );

  throwExceptionOnFailure( checkIfEqual( s11, 0.5 * 1 / 2 * 200 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
  throwExceptionOnFailure( checkIfEqual( s22, 1. / 2 * 200 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
  throwExceptionOnFailure( checkIfEqual( s33, 0.2 * 1 / 2 * 200 ), MakeString() << __PRETTY_FUNCTION__ << " failed" );
}

int main()
{

  test_geostaticStress();

  return 0;
}
