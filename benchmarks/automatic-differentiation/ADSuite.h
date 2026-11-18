#pragma once
#include <functional>
#include <string>
#include <vector>

namespace ADBench {

  struct Case {
    std::string                              id;       // kernel family; same id = same outputs across suites
    std::string                              name;     // row label
    std::function< std::vector< double >() > evaluate; // one kernel call, all outputs flattened
  };

  struct Suite {
    std::string         label; // backend + mode, e.g. "sacado-vec"
    std::vector< Case > cases;
  };

  std::vector< Suite > makeAutodiffSuites();
  std::vector< Suite > makeBoostSuites();
  std::vector< Suite > makeSacadoSuites();
  std::vector< Suite > makeCeresSuites();
  std::vector< Suite > makeStanSuites();
  std::vector< Suite > makeTinyADSuites();
  std::vector< Suite > makeAdeptSuites();
  std::vector< Suite > makeMarmotSuites();

} // namespace ADBench
