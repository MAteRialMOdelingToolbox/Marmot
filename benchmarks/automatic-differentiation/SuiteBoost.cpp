#include "ADBackendBoost.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeBoostSuites()
{
  return makeSuites< BoostBackend >();
}
