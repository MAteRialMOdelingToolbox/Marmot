#include "Marmot/MarmotAutomaticDifferentiation.h"
#include <autodiff/forward/utils/derivative.hpp>
#include <iostream>

using namespace autodiff;
using namespace Eigen;

namespace Marmot {

  namespace AutomaticDifferentiation {

    dual2nd shiftTo2ndOrderDual( const dual& x )
    {
      dual2nd x2nd( 0.0 );
      x2nd.val = x;

      return x2nd;
    }

    double df_dx( const scalar_to_scalar_function_type& f, const double& x )
    {
      dual x_right;
      x_right.val = x;
      seed< 1 >( x_right, 1.0 );

      const double df_dx = f( x_right ).grad;

      return df_dx;
    }

    dual df_dx( const scalar_to_scalar_function_type_2nd& f, const dual& x )
    {

      dual2nd x_right = shiftTo2ndOrderDual( x );
      seed< 1 >( x_right, 1.0 );

      dual          df_dx;
      const dual2nd f_right = f( x_right );
      df_dx.val             = derivative< 1 >( f_right );
      df_dx.grad            = derivative< 2 >( f_right );

      return df_dx;
    }
    MatrixXd forwardMode( const vector_to_vector_function_type& F, const VectorXd& X )
    {
      VectorXdual X_( X );
      const auto  J = jacobian( F, wrt( X_ ), at( X_ ) );

      return J;
    }

    std::pair< VectorXd, MatrixXd > jacobian( const vector_to_vector_function_type_dual& F, const VectorXd& X )
    {
      VectorXdual  X_( X );
      VectorXdual  F_right;
      const size_t sizeX = X_.rows();
      MatrixXd     J( sizeX, sizeX );
      VectorXd     F_( sizeX );

      // J_ij = d F_i / d x_j

      for ( size_t j = 0; j < sizeX; j++ ) {

        seed< 1 >( X_( j ), 1.0 );
        F_right = F( X_ );

        for ( size_t i = 0; i < sizeX; i++ ) {
          J( i, j ) = derivative< 1 >( F_right( i ) );
        }
        seed< 1 >( X_( j ), 0.0 );
        F_( j ) = F_right( j ).val;
      }

      return { F_, J };
    }
  } // namespace AutomaticDifferentiation
} // namespace Marmot
