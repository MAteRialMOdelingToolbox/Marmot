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

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > dF_dX( const function_type& F, const Eigen::VectorXd& X )
    {
      using namespace autodiff;
      VectorXdual     X_( X );
      const size_t    rows = X.rows();
      VectorXdual     FDual( rows );
      Eigen::VectorXd F_( rows );

      const Eigen::MatrixXd J = jacobian( F, wrt( X_ ), at( X_ ), FDual );

      for ( size_t i = 0; i < rows; i++ ) {
        F_( i ) = FDual( i ).val;
      }

      return { F_, J };
    }
  } // namespace AutomaticDifferentiation
} // namespace Marmot
