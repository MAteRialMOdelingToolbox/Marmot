#include "Marmot/MarmotTesting.h"
#include <iostream>

namespace Marmot::Testing {

  std::string printScalar( const double a )
  {
    return std::to_string( a );
  }

  std::string printScalar( const autodiff::dual a )
  {
    return "(" + std::to_string( a.val ) + ",  " + std::to_string( a.grad ) + ")";
  }

  bool checkIfEqual( const double a, const double b, const double tol )
  {
    if ( std::abs( a - b ) > tol ) {
      std::cout << " Hint: a = " << a << " != "
                << " b = " << b << " ( tol = " << tol << " )" << std::endl;
      return false;
    }
    else
      return true;
  }

  bool checkIfEqual( const autodiff::dual a, const autodiff::dual b, const double tol )
  {
    if ( checkIfEqual( a.val, b.val, tol ) && checkIfEqual( a.grad, b.grad, tol ) ) {
      return true;
    }
    return false;
  }

  void throwExceptionOnFailure( const bool condition, const std::string& message )
  {
    if ( !condition ) {
      throw std::runtime_error( message );
    }
  }

} // namespace Marmot::Testing
