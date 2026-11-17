#include "Marmot/MarmotAutomaticDifferentiation.h"

namespace Marmot::Math::Differentiation::Automatic {

  dual2nd shiftTo2ndOrderDual( const dual& x )
  {
    return increaseDualOrderWithShift< 1 >( x );
  }

  VectorXdual2nd shiftTo2ndOrderDual( const VectorXdual& X )
  {
    const size_t   sizeX = X.size();
    VectorXdual2nd X_( sizeX );

    for ( size_t j = 0; j < sizeX; j++ ) {
      X_( j ) = shiftTo2ndOrderDual( X( j ) );
    }

    return X_;
  }

  double df_dx( const scalar_to_scalar_function_type& f, const double& x )
  {
    dual x_right( x );
    seed< 1 >( x_right, 1.0 );

    const double df_dx = derivative< 1 >( f( x_right ) );

    return df_dx;
  }

  dual df_dx( const scalar_to_scalar_function_type_2nd& f, const dual& x )
  {

    dual2nd x_right = shiftTo2ndOrderDual( x );
    seed< 1 >( x_right, 1.0 );

    dual2nd f_right = f( x_right );

    // ( derivative< 1 >, derivative< 2 > ) of f_right
    return decreaseDualOrderWithShift< 2 >( f_right );
  }

  std::pair< Eigen::VectorXd, Eigen::MatrixXd > dF_dX( const vector_to_vector_function_type_dual& F,
                                                       const Eigen::VectorXd&                     X )
  {
    VectorXdual     X_ = X.cast< dual >();
    VectorXdual     F_right;
    const size_t    sizeX = X_.rows();
    Eigen::MatrixXd J( sizeX, sizeX );
    Eigen::VectorXd F_( sizeX );

    // J_ij = d F_i / d x_j
    for ( size_t j = 0; j < sizeX; j++ ) {

      seed< 1 >( X_( j ), 1.0 );
      F_right = F( X_ );

      for ( size_t i = 0; i < sizeX; i++ ) {
        J( i, j ) = derivative< 1 >( F_right( i ) );
      }
      seed< 1 >( X_( j ), 0.0 );
      F_( j ) = value( F_right( j ) );
    }

    return { F_, J };
  }

  std::pair< VectorXdual, MatrixXdual > dF_dX_2nd( const vector_to_vector_function_type_dual2nd& F,
                                                   const VectorXdual&                            X )
  {
    VectorXdual2nd X_ = increaseDualOrderWithShift< 1 >( X );
    VectorXdual2nd F_right;
    const size_t   sizeX = X_.rows();
    MatrixXdual    J( sizeX, sizeX );
    VectorXdual    F_( sizeX );

    // J_ij = d F_i / d x_j
    for ( size_t j = 0; j < sizeX; j++ ) {

      seed< 1 >( X_( j ), 1.0 );
      F_right = F( X_ );

      for ( size_t i = 0; i < sizeX; i++ ) {
        // ( derivative< 1 >, derivative< 2 > ) of F_right( i )
        J( i, j ) = decreaseDualOrderWithShift< 2 >( F_right( i ) );
      }

      seed< 1 >( X_( j ), 0.0 );
      F_( j ) = dual( value( F_right( j ) ) );
      seed< 1 >( F_( j ), derivative< 1 >( F_right( j ) ) );
    }

    return { F_, J };
  }
} // namespace Marmot::Math::Differentiation::Automatic
