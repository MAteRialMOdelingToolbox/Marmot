#include "ADBackendAutodiff.h"
#include "ADSuiteImpl.h"

ADBench::Suite ADBench::makeAutodiffSuite()
{
  return makeSuite< AutodiffBackend >();
}
