#include "ADBackendStan.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeStanSuites()
{
  return makeSuites< StanBackend >();
}
