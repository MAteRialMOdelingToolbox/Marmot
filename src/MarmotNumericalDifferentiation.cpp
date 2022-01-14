#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include <complex>

using namespace Eigen;

namespace Marmot {
  namespace NumericalAlgorithms::Differentiation {

    MatrixXd forwardDifference( const function_type& F, const VectorXd& X )
    {

      const auto xSize = X.rows();
      MatrixXd   J( xSize, xSize );

      VectorXd rightX( xSize );

      for ( auto i = 0; i < xSize; i++ ) {
        double volatile h = std::max( 1.0, std::abs( X( i ) ) ) * Marmot::Constants::SquareRootEps;
        // clang-format off
        rightX = X;
        rightX( i ) += h;

        J.col( i ) = (  F( rightX )  - F( X  ) )
            / //------------------------------------
                            ( 1. * h );
        // clang-format on
      }

      return J;
    }

    MatrixXd centralDifference( const function_type& F, const VectorXd& X )
    {

      const auto xSize = X.rows();
      MatrixXd   J( xSize, xSize );

      VectorXd leftX( xSize );
      VectorXd rightX( xSize );

      for ( auto i = 0; i < xSize; i++ ) {
        double volatile h = std::max( 1.0, std::abs( X( i ) ) ) * Marmot::Constants::CubicRootEps;
        // clang-format off
        leftX  = X;
        rightX = X;
        leftX( i ) -= h;
        rightX( i ) += h;

        J.col( i ) = (  F( rightX )  - F( leftX ) )
            / //------------------------------------
                             ( 2. * h );
        // clang-format on
      }

      return J;
    }

    namespace Complex {
      /*
       * Implementation of Numerical Differantiation using Complex Step Approximations
       *
       * further Information can be found in
       *   - Martins et al. (2003) The Complex-Step Derivative Approximation
       *   - Lai et al. (2005) New Complex-Step Derivative Approximations ...
       *
       */

      const static std::complex< double > imaginaryUnit = { 0, 1 };
      const static std::complex< double > complexUnit   = { 1, 1 };
      const static std::complex< double > i_            = Marmot::Constants::sqrt2 / 2. * complexUnit;

      MatrixXd forwardDifference( const function_type_complex& F, const VectorXd& X )
      {
        /*
         * according to Martins et al. (2003) Equ. 6
         * according to Lai et al. (2005) Equ. 7
         */
        const auto xSize = X.rows();
        MatrixXd   J( xSize, xSize );
        VectorXcd  rightX( xSize );

        for ( auto i = 0; i < xSize; i++ ) {
          double h = std::max( 1.0, std::abs( X( i ) ) ) * 1e-16;
          rightX   = X;
          rightX( i ) += h * imaginaryUnit;

          J.col( i ) = F( rightX ).imag() / h;
        }

        return J;
      }

      MatrixXd centralDifference( const function_type_complex& F, const VectorXd& X )
      {
        /*
         * according to Lai et al. (2005) Equ. 19
         */
        const auto xSize = X.rows();
        MatrixXd   J( xSize, xSize );
        VectorXcd  e( xSize );
        VectorXcd  rightX( xSize );
        VectorXcd  leftX( xSize );

        complexDouble i_ = Marmot::Constants::sqrt2 / 2. * complexUnit;

        for ( auto i = 0; i < xSize; i++ ) {
          double h = std::max( 1.0, std::abs( X( i ) ) ) * Marmot::Constants::SquareRootEps;

          e.setZero();
          e( i ) = 1.;

          leftX = X;
          leftX -= e * i_ * h;

          rightX = X;
          rightX += e * i_ * h;

          // clang-format off
          J.col( i ) =      ( F( rightX ) - F( leftX )  ).imag() 
                       / //--------------------------------------
                              ( Marmot::Constants::sqrt2 * h );

          // clang-format on
        }
        return J;
      }

      MatrixXd fourthOrderAccurateDerivative( const function_type_complex& F, const VectorXd& X )
      {
        /*
         * according to Lai et al. (2005) Equ. 24
         */
        const auto xSize = X.rows();
        MatrixXd   J( xSize, xSize );
        VectorXcd  complex_X( X );
        VectorXcd  e( xSize );
        VectorXcd  x1_( xSize );
        VectorXcd  x2_( xSize );
        VectorXcd  x3_( xSize );
        VectorXcd  x4_( xSize );

        for ( auto i = 0; i < xSize; i++ ) {
          double h = std::max( 1.0, std::abs( X( i ) ) ) * Marmot::Constants::SquareRootEps;
          e.setZero();
          e( i ) = 1.;

          x1_ = complex_X + e / 2. * i_ * h;
          x2_ = complex_X - e / 2. * i_ * h;
          x3_ = complex_X + e * i_ * h;
          x4_ = complex_X - e * i_ * h;

          // clang-format off
          J.col( i ) = ( 8.* ( F( x1_ ) - F( x2_ ) ) 
                           - ( F( x3_ ) - F( x4_ ) ) ).imag() 
                           / ( Marmot::Constants::sqrt2 * 3. * h );

          // clang-format on
        }

        return J;
      }
    } // namespace Complex

  } // namespace NumericalAlgorithms::Differentiation
} // namespace Marmot
