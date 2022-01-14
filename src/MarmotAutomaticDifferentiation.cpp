#include "Marmot/MarmotAutomaticDifferentiation.h"
#include <iostream>

namespace Marmot {
  namespace AutomaticDifferentiation {

    Eigen::MatrixXd forwardMode( const function_type& F, const Eigen::VectorXd& X )
    {
      using namespace autodiff;
      VectorXdual X_( X );
      const auto  J = jacobian( F, wrt( X_ ), at( X_ ) );

      return J;
    }

  } // namespace AutomaticDifferentiation
} // namespace Marmot
