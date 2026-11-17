#include "ADBackendTinyAD.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeTinyADSuites()
{
  return makeSuites< TinyADBackend >();
}
