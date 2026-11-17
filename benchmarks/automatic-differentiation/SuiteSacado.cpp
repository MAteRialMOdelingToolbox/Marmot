#include "ADBackendSacado.h"
#include "ADSuiteImpl.h"

std::vector< ADBench::Suite > ADBench::makeSacadoSuites()
{
  return makeSuites< SacadoBackend >();
}
