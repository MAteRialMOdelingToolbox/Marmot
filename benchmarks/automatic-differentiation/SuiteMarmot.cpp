#include "ADBackendMarmot.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeMarmotSuites()
{
  return makeSuites< MarmotBackend >();
}
