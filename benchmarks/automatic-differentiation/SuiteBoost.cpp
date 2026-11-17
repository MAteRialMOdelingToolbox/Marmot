#include "ADBackendBoost.h"
#include "ADSuiteImpl.h"

ADBench::Suite ADBench::makeBoostSuite()
{
  return makeSuite< BoostBackend >();
}
