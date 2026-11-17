#pragma once
#include <functional>
#include <string>
#include <vector>

namespace ADBench {

  struct Case {
    std::string                             name;
    std::string                             mirrors;  // Marmot code path the kernel stands for
    std::function< std::vector< double >() > evaluate; // one kernel call, all outputs flattened
  };

  struct Suite {
    std::string         backend;
    std::vector< Case > cases;
  };

  Suite makeAutodiffSuite();
  Suite makeBoostSuite();

} // namespace ADBench
