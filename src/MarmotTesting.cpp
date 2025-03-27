#include "Marmot/MarmotTesting.h"
#include <iostream>

namespace Marmot::Testing {

  std::string getString( const double a )
  {
    return std::to_string( a );
  }

  std::string getString( const autodiff::dual a )
  {
    return "(" + std::to_string( a.val ) + ",  " + std::to_string( a.grad ) + ")";
  }

  bool checkIfEqual( const double a, const double b, const double tol )
  {
    // check for NaN values and return false if any of the values is NaN
    if ( std::isnan( a ) || std::isnan( b ) ) {
      std::cout << " Hint: a = " << a << " or "
                << " b = " << b << " is NaN" << std::endl;
      return false;
    }

    // check for inf values and return false if any of the values is inf
    if ( std::isinf( a ) || std::isinf( b ) ) {
      std::cout << " Hint: a = " << a << " or "
                << " b = " << b << " is inf" << std::endl;
      return false;
    }

    // check if values are equal with respect to the tolerance
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
