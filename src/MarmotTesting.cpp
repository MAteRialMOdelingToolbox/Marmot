#include <exception>
#include <iostream>
#include <vector>

#include "Marmot/MarmotTesting.h"

namespace Marmot::Testing {

  void checkIfEqual( double a, double b )
  {

    if ( std::abs( a - b ) > 1e-15 ) {
      throw std::exception();
    };
  }

  void checkIfEqual( autodiff::dual a, autodiff::dual b )
  {

    if ( std::abs( a.val - b.val ) > 1e-15 || std::abs( a.grad - b.grad ) > 1e-15 ) {
      throw std::exception();
    };
  }
} // namespace Marmot::Testing
