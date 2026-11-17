#include "ADBackendAutodiff.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeAutodiffSuites()
{
  return makeSuites< AutodiffBackend >();
}
