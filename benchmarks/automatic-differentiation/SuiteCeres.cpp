#include "ADBackendCeres.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeCeresSuites()
{
  return makeSuites< CeresBackend >();
}
