#include "ADBackendAdept.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeAdeptSuites()
{
  return makeSuites< AdeptBackend >();
}
